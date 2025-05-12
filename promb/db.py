from __future__ import annotations
import gzip
import os
import sys
import time
from typing import Sequence, Iterable

from promb.__version__ import __version__
import heapq
import itertools
import math
from tqdm import tqdm
from promb.parsing import parse_fasta
import pandas as pd
import numpy as np


PSEUDOPROB_EPSILON = 1e-30
HUMAN_OASIS_DB_PATH = os.path.join(os.path.dirname(__file__), 'resources', 'OASis_9mers_v1_10perc_subjects.txt.gz')
HUMAN_REFERENCE_DB_PATH = os.path.join(os.path.dirname(__file__), 'resources', 'uniprotkb_proteome_UP000005640_2025_04_21.fasta.gz')
HUMAN_SWISSPROT_DB_PATH = os.path.join(os.path.dirname(__file__), 'resources', 'uniprot_sprot_human.fasta.gz')

TEXT_EXTENSIONS = ('.txt', '.txt.gz')
FASTA_EXTENSIONS = ('.fa.gz', '.fasta.gz', '.fa', '.fasta')
DB_EXTENSIONS = TEXT_EXTENSIONS + FASTA_EXTENSIONS

BLOSUM62 = pd.DataFrame([[4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0, 0.0, -3.0, -2.0, 0.0, -2.0, -1.0, 0.0, -4.0], [-1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0, -1.0, -1.0, -3.0, -2.0, -3.0, -1.0, 0.0, -1.0, -4.0], [-2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0, 1.0, 0.0, -4.0, -2.0, -3.0, 3.0, 0.0, -1.0, -4.0], [-2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0, 0.0, -1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0], [0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -1.0, -2.0, -2.0, -1.0, -3.0, -3.0, -2.0, -4.0], [-1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0, -1.0, -2.0, -1.0, -2.0, 0.0, 3.0, -1.0, -4.0], [-1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0], [0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0, 0.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -1.0, -4.0], [-2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0, -2.0, -2.0, 2.0, -3.0, 0.0, 0.0, -1.0, -4.0], [-1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0, -2.0, -1.0, -3.0, -1.0, 3.0, -3.0, -3.0, -1.0, -4.0], [-1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0, -2.0, -1.0, -2.0, -1.0, 1.0, -4.0, -3.0, -1.0, -4.0], [-1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0, 0.0, 1.0, -1.0, -4.0], [-1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, 1.0, -3.0, -1.0, -1.0, -4.0], [-2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0, -2.0, -2.0, 1.0, 3.0, -1.0, -3.0, -3.0, -1.0, -4.0], [-1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 7.0, -1.0, -1.0, -4.0, -3.0, -2.0, -2.0, -1.0, -2.0, -4.0], [1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0, 4.0, 1.0, -3.0, -2.0, -2.0, 0.0, 0.0, 0.0, -4.0], [0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0, 5.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -4.0], [-3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0, -4.0, -3.0, -2.0, 11.0, 2.0, -3.0, -4.0, -3.0, -2.0, -4.0], [-2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0, -3.0, -2.0, -2.0, 2.0, 7.0, -1.0, -3.0, -2.0, -1.0, -4.0], [0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0, -2.0, 0.0, -3.0, -1.0, 4.0, -3.0, -2.0, -1.0, -4.0], [-2.0, -1.0, 3.0, 4.0, -3.0, 0.0, 1.0, -1.0, 0.0, -3.0, -4.0, 0.0, -3.0, -3.0, -2.0, 0.0, -1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0], [-1.0, 0.0, 0.0, 1.0, -3.0, 3.0, 4.0, -2.0, 0.0, -3.0, -3.0, 1.0, -1.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0], [0.0, -1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0, 0.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -4.0], [-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0]])
BLOSUM62.index = list('ARNDCQEGHILKMFPSTWYVBZX*')
BLOSUM62.columns = list('ARNDCQEGHILKMFPSTWYVBZX*')
BLOSUM62 = BLOSUM62.to_dict(orient='index')

def init_db(db_path: str, peptide_length=None, ignored_peptides: Sequence[str] = None, ignored_ids: Sequence[str] = None, verbose=True) -> 'PrombDB':
    """Initialize database from given path (.txt.gz, one peptide per line)

    Pre-packaged databases are available as constants: "human-oas", "human-swissprot"

    peptide_length needs to be provided when loading a .fasta.gz database
    """
    if db_path is None:
        raise ValueError('Use init_db("human-oas") or init_db("human-swissprot", 9) or init_db("path/to/database.fasta.gz", 9)')
    if db_path.lower() == "human-oas":
        db_path = HUMAN_OASIS_DB_PATH
    elif db_path.lower() == "human-reference":
        db_path = HUMAN_REFERENCE_DB_PATH
    elif db_path.lower() == "human-swissprot":
        db_path = HUMAN_SWISSPROT_DB_PATH
    elif '.' not in db_path:
        raise ValueError(f'Unknown database {db_path}, use "human-reference", "human-swissprot", "human-oas", or provide a file path (txt or fasta)')
    # Load database from disk
    if db_path.endswith(TEXT_EXTENSIONS):
        if peptide_length is not None:
            raise ValueError('txt databases should contain peptides of the same length (one peptide per line), peptide_length is not needed')
    elif db_path.endswith(FASTA_EXTENSIONS):
        if peptide_length is None:
            raise ValueError('Please provide peptide_length when reading fasta databases')
    else:
        raise ValueError(f'Unsupported DB file format: {db_path}, supported formats: {DB_EXTENSIONS}')
    if verbose:
        print('''
              *           █
▄▄▄▄    ▄▄▄  ▄█▄   ▄▄▄▄   █▄▄▄  
█   █  █    █▓███  █ █ █  █   █ 
█▄▄▄▀  █    ▀███▀  █   █  █▄▄▄▀ 
█                              
▀  protein mutation burden v{}
                                '''.format(__version__), file=sys.stderr)
        print('Loading promb database to memory...', file=sys.stderr)
        start_time = time.time()

    with (gzip.open(db_path, 'rt') if db_path.endswith('.gz') else open(db_path, 'rt')) as f:
        if db_path.endswith(TEXT_EXTENSIONS):
            # read peptides from file
            peptides = frozenset(line.strip() for line in f)
            sequences = None
        elif db_path.endswith(FASTA_EXTENSIONS):
            # read peptides from fasta file
            sequences = parse_fasta(f)
            if ignored_ids is not None:
                for ignored_id in ignored_ids:
                    del sequences[ignored_id] 
            peptides = frozenset(seq[i : i + peptide_length] for seq in sequences.values() for i in range(len(seq) - peptide_length + 1))
        else:
            raise ValueError(f'Unsupported DB file format: {db_path}, supported formats: {DB_EXTENSIONS}')

        if ignored_peptides:
            peptides = frozenset(peptides - set(ignored_peptides))

        if verbose:
            print(f'Database '
                  f'(~{len(peptides) / 1000000:.0f}M peptides, <{sys.getsizeof(peptides) / 1024 / 1024:.0f} MB) '
                  f'loaded in {time.time() - start_time:.0f}s.', file=sys.stderr)
    return PrombDB(peptides, sequences=sequences)


class PrombDB:
    def __init__(self, peptides: frozenset[str], sequences: dict[str, str] | None = None):
        assert isinstance(peptides, frozenset), f'Expected frozenset, got {type(peptides)}'
        peptide_lengths = sorted(set(len(peptide) for peptide in peptides))
        assert len(peptide_lengths) == 1, f'All peptides in a DB should be the same length, found lengths: {peptide_lengths}'
        self.peptide_length = peptide_lengths[0]
        self.peptides = peptides
        self.sequences = sequences
        self._db_peptides_encoded = None

    def __contains__(self, peptide: str) -> bool:
        return self.contains(peptide)

    def contains(self, peptide: str) -> bool:
        """Check if peptide is present (exact match)."""
        assert len(peptide) == self.peptide_length, f"Expected peptide of length {self.peptide_length}, got sequence: {peptide}"
        return peptide in self.peptides

    def contains_or_point_mutant(self, peptide: str) -> bool:
        """Check if peptide is in reference or has at most one mutation from the closest reference peptide"""
        return self.contains(peptide) or any(self.find_point_mutant_peptides_iter(peptide))

    def contains_or_double_mutant(self, peptide: str) -> bool:
        """Check if peptide is in reference or has at most two mutations from a the closest reference peptide"""
        return self.contains_or_point_mutant(peptide) or any(self.find_double_mutant_peptides_iter(peptide))
    
    def find_peptide_sources(self, peptide: str) -> dict[str, str]:
        """Find sequences in our DB that contain a given peptide, return dict (header -> seq)"""
        assert self.sequences is not None, 'Finding peptide source is only supported with full sequence-based DBs'
        return {header: seq for header, seq in self.sequences.items() if peptide in seq}

    def compute_peptide_content(self, seq: str) -> float:
        """Compute peptide content, for example human peptide content: Fraction of peptides that are completely human (exact match)"""
        query_peptides = self.chop_seq_peptides(seq)
        num_human_peptides = sum(self.contains(peptide) for peptide in query_peptides)
        return num_human_peptides / len(query_peptides)

    def compute_pmb(self, seq: str, max: int = None) -> float:
        """Compute PMB: Protein mutation burden (average number of mutations to closest peptide in reference DB)
        
        When max is provided, give up searching at this number of mutations and return max for the given peptide.
        """
        query_peptides = self.chop_seq_peptides(seq)
        peptide_wise_mutations = self.compute_peptide_wise_mutations(query_peptides, max=max)
        return sum(peptide_wise_mutations) / len(query_peptides)

    def compute_peptide_wise_mutations(self, query_peptides: list[str], max: int = None) -> list[float]:
        """Compute peptide-wise PMB: Peptide mutation burden (number of mutations to closest peptide in reference DB)
        
        When max is provided, give up searching at this number of mutations and return max for the given peptide.
        """
        assert not isinstance(query_peptides, str), "compute_peptide_wise_mutations expects a list of peptides, please use db.chop_seq_peptides(seq)"
        if not query_peptides:
            return []
        if max:
            assert 1 <= max <= 3, 'Max should be None, 1, 2, or 3'
        results = []
        multiple_mutants = []
        multiple_mutant_indexes = []
        for i, peptide in enumerate(query_peptides):
            assert isinstance(peptide, str) and len(peptide) == self.peptide_length, \
                f'Expected list of {self.peptide_length}-mer strings, got: {peptide}, ...'
            # Speed up scoring by finding human peptides, point mutants and double mutants (quick)
            if self.contains(peptide):
                results.append(0)
            elif (max and max <= 1) or self.contains_or_point_mutant(peptide):
                results.append(1)
            elif (max and max <= 2) or self.contains_or_double_mutant(peptide):
                results.append(2)
            elif max:
                if max <= 3:
                    results.append(3)
                else:
                    raise ValueError()
            else:
                multiple_mutants.append(peptide)
                multiple_mutant_indexes.append(i)
                results.append(None)
        nearest_peptides = self.find_nearest_peptides(multiple_mutants, 1)
        for i, peptide, nearest in zip(multiple_mutant_indexes, multiple_mutants, nearest_peptides):
            results[i] = sum(aa != bb for aa, bb in zip(peptide, nearest[0]))
        return results

    def chop_seq_peptides(self, seq) -> list[str]:
        """Chop sequence into overlapping peptides of given length"""
        assert isinstance(seq, str), f'Expected string, got: {type(seq)} {seq}'
        return [seq[i : i + self.peptide_length] for i in range(len(seq) - self.peptide_length + 1)]

    def find_nearest_peptides(self, query_peptides: list[str], n: int) -> list[list[str]]:
        """Find nearest peptides in the reference DB for a list of peptides, returns list of lists
        
        Peptides are ordered by number of mutations, and then (in case of ties) by BLOSUM similarity score.
        """
        try:
            from scipy.spatial.distance import cdist
            import numpy as np
        except ImportError:
            raise ValueError('Please install "scipy" module to use this function.')
        
        assert not isinstance(query_peptides, str), "find_nearest_peptides expects a list of peptides, please use db.chop_seq_peptides(seq)"
    
        if not query_peptides:
            return []

        # we need to sort the list to make sure we get the same result each time
        all_peptides = np.array(sorted(self.peptides))
        if self._db_peptides_encoded is None:
            self._db_peptides_encoded = np.array([[ord(aa) for aa in peptide] for peptide in all_peptides], dtype='double')
        assert len(self._db_peptides_encoded) == len(self.peptides)

        nearest_peptides = []
        for peptide in tqdm(query_peptides):
            if n == 1 and self.contains(peptide):
                nearest_peptides.append([peptide])
                continue
            distances = cdist([[ord(aa) for aa in peptide]], self._db_peptides_encoded, metric='hamming')[0] * self.peptide_length
            # find peptides with up to 0, 1, 2, ... mutations until we get at least n hits
            for num_mutations in range(0, self.peptide_length):
                is_hit = (distances <= num_mutations + 0.000001)
                if is_hit.sum() >= n:
                    break
            hits = all_peptides[is_hit]
            if len(hits) == 1:
                # we don't need to sort
                nearest_peptides.append([str(hits[0])])
                continue
            # sort hits by identity and then resolve ties by BLOSUM distance
            distances2 = [distance * 10000 - sum(BLOSUM62.get(aa, {}).get(bb, -10) for aa, bb in zip(peptide, hit)) for hit, distance in zip(hits, distances[is_hit])]
            hits = hits[np.argsort(distances2, kind='stable')[:n]]
            nearest_peptides.append([str(p) for p in hits])
        return nearest_peptides

    def find_point_mutant_peptides(self, peptide: str, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> list[str]:
        """Find reference peptides with at most one mutation from provided sequence"""
        return list(self.find_point_mutant_peptides_iter(peptide, allowed_amino_acids=allowed_amino_acids))
    
    def find_point_mutant_peptides_iter(self, peptide: str, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> Iterable[str]:
        return (mutant
                for offset in range(self.peptide_length)
                for mutant in self.find_peptides_matching_template(add_wildcard(peptide, offset), allowed_amino_acids=allowed_amino_acids)
                if mutant != peptide)
    
    def find_double_mutant_peptides(self, peptide: str, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> list[str]:
        """Find reference peptides with at most two mutations from provided sequence"""
        return list(self.find_double_mutant_peptides_iter(peptide, allowed_amino_acids=allowed_amino_acids))
    
    def find_double_mutant_peptides_iter(self, peptide: str, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> Iterable[str]:
        return (mutant
                for offset1 in range(self.peptide_length)
                for offset2 in range(self.peptide_length)
                for mutant in self.find_peptides_matching_template(add_wildcard(add_wildcard(peptide, offset1), offset2), allowed_amino_acids=allowed_amino_acids)
                if mutant != peptide and offset1 != offset2)
    
    def find_peptides_matching_template(self, peptide: str, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> list[str]:
        """Find reference peptides matching sequence and wildcards (*)
    
        NOTE: This function increases in complexity with every wildcard.
    
        Example: find_peptides_matching_template('QVQL*QSGA') -> ['QVQLVQSGA', 'QVQLLQSGA', ...]
        """
        if '*' not in peptide:
            # optimization, check for exact match
            return [peptide] if self.contains(peptide) else []
        assert len(peptide) == self.peptide_length, f"Expected peptide of length {self.peptide_length}, got sequence: {peptide}"
        if peptide.count('*') == 1:
            # 20 combinations, naive implementation
            peptides = []
            for aa in allowed_amino_acids:
                new_peptide = peptide.replace('*', aa)
                if self.contains(new_peptide):
                    peptides.append(new_peptide)
            return peptides
        aa_combinations = itertools.product(allowed_amino_acids, repeat=peptide.count('*'))
        peptides = []
        for aa_combination in aa_combinations:
            new_peptide = peptide
            for aa in aa_combination:
                new_peptide = new_peptide.replace('*', aa, 1)
            if self.contains(new_peptide):
                peptides.append(new_peptide)
        return peptides
    
    def find_peptides_matching_distribution(self, distribution: list[list[float]], alphabet: str | list[str] | dict[str, int], n: int) -> list[str]:
        """Find reference peptides with highest probability according to distribution
    
        :param distribution: Table of probabilities (9 positions x 20 amino acids - or depending on alphabet)
        :param alphabet: List of amino acids
        :param n: Number of peptides to return
    
        :return: List of n peptides with highest probability
        """
    
        try:
            first_row = distribution[0]
        except:
            raise ValueError(f'Expected matrix 9 x {len(alphabet)}, got flat array')
    
        assert len(distribution) == self.peptide_length, f'Expected matrix 9 x {len(alphabet)} (position x amino acid), got {len(distribution)} rows'
        assert len(distribution[0]) == len(alphabet), f'Expected {len(alphabet)} probabilities, one for each amino acid, got: {len(distribution)}'
        if not isinstance(alphabet, dict):
            alphabet = {aa: i for i, aa in enumerate(alphabet)}
        log_distribution = [[math.log10(prob + PSEUDOPROB_EPSILON) for prob in probs] for probs in distribution]
        return heapq.nlargest(
            n=n,
            iterable=self.peptides,
            key=lambda peptide: sum(log_probs[alphabet[aa]] if aa in alphabet else math.log10(PSEUDOPROB_EPSILON) for log_probs, aa in zip(log_distribution, peptide))
        )

    def compute_pssm(self, seq: str, nearest_peptides=1, normalize=True, ignore_wildtype=False, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> pd.DataFrame:
        """Compute position-specific scoring matrix (PSSM) for a sequence - fraction of nearest peptides in the reference database
        containing a given amino acid (column) at a given position (row)

        NOTE this does NOT take into account the frequency of the hit peptides (how many sequences they appear in)
        but rather the number of nearest peptides.

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param normalize: Normalize to frequencies at each position. When false, returns counts (of nearest peptides containing a given mutation).
        :param ignore_wildtype: Do not count amino acids present in the input sequence at the given position.
        :return: DataFrame with mutation counts at each position (rows) and amino acids (columns)
        """
        peptides = self.chop_seq_peptides(seq)
        # pos -> aa -> count
        counts = {}
        if isinstance(nearest_peptides, list):
            assert len(nearest_peptides) == len(peptides), f"Expected list of length {len(peptides)} for sequence of length {len(seq)}, got list of length {len(nearest_peptides)}"
        elif isinstance(nearest_peptides, int):
            nearest_peptides = self.find_nearest_peptides(peptides, n=nearest_peptides)
        else:
            raise ValueError("nearest_peptides should be int (number of nearest peptides) or list of lists (nearest peptide hits)")
        for start, (query, hits) in enumerate(zip(peptides, nearest_peptides)):
            for hit in hits:
                for pos, aa in enumerate(hit):
                    if start + pos + 1 not in counts:
                        counts[start + pos + 1] = {}
                    counts[start + pos + 1][aa] = counts[start + pos + 1].get(aa, 0) + 1
        counts = pd.DataFrame(counts).reindex(list(allowed_amino_acids)).T.reindex(range(1, len(seq) + 1)).fillna(0)
        counts.index.name = 'position'
        total_counts_by_pos = counts.sum(axis=1)
        if ignore_wildtype:
            for i, aa in enumerate(seq, start=1):
                counts.loc[i, aa] = 0
        if not normalize:
            return counts
        return counts.div(total_counts_by_pos, axis=0).fillna(0)

    def compute_positional_likelihood(self, seq: str, nearest_peptides=1) -> list[float]:
        """Compute Fraction of nearest overlapping {peptide_length}mers that contain the input amino acid at that position.

        Positions with values close to 0 can be considered non-human, and values close to 1 can be considered human

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find (must be 1), or list of lists of nearest peptides from db.find_nearest_peptides()
        :return: list of likelihoods, one for each positions
        """
        if isinstance(nearest_peptides, list):
            if len(nearest_peptides[0]) != 1:
                print("Note: Using only nearest 1 peptide to compute positional likelihood to ensure amino acids in peptides already found in the DB (already human) do not get flagged", file=sys.stderr)
                nearest_peptides = [hits[:1] for hits in nearest_peptides]
        elif isinstance(nearest_peptides, int):
            assert nearest_peptides == 1, "compute_positional_likelihood should use 1 nearest peptide (to ensure peptides already found in the DB do not get flagged)"
        else:
            raise ValueError("nearest_peptides should be int (number of nearest peptides) or list of lists (nearest peptide hits)")
        mutation_freqs = self.compute_pssm(
            seq,
            nearest_peptides=nearest_peptides,
            ignore_wildtype=True
        )
        return (1 - mutation_freqs.sum(axis=1)).tolist()

    def suggest_point_mutant_candidates(self, seq: str, nearest_peptides=1, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> list[str]:
        """Suggest point mutant sequences for a given sequence to reduce PMB (avg distance to nearest reference peptides)

        NOTE that this does not guarantee complete humanness - but it should
        reduce the "PMB" - average distance to nearest human peptide.

        ALSO NOTE this does NOT take into account the frequency of the hit peptides (how many sequences they appear in)
        but rather the number of unique peptides in the reference db that contain the given point mutation.

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param allowed_amino_acids: Allowed amino acids for mutations
        :return: List of point mutant sequences, from most frequent
        """
        pssm_mutations = self.compute_pssm(seq, nearest_peptides=nearest_peptides, ignore_wildtype=True, allowed_amino_acids=allowed_amino_acids)
        sorted_freqs = pssm_mutations.unstack().sort_values(ascending=False)

        # Iterate and print
        mutants = []
        for (mut_aa, pos), val in sorted_freqs.items():
            if val == 0:
                break
            mutant = list(seq)
            mutant[pos-1] = mut_aa
            mutants.append(''.join(mutant))
        return mutants

    def find_point_mutant_seqs(self, seq: str, nearest_peptides=1, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY'):
        """Exhaustively find all possible point mutant sequences for a given sequence where all overlapping peptides
        are found in the reference database (are human).

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param allowed_amino_acids: Allowed amino acids for mutations
        :return: List of point mutant sequences (if any) whose peptides are contained in the db
        """
        candidates = self.suggest_point_mutant_candidates(seq, nearest_peptides=nearest_peptides, allowed_amino_acids=allowed_amino_acids)
        return [candidate for candidate in candidates if self.compute_peptide_content(candidate) == 1.0]

    def suggest_double_mutant_candidates(self, seq: str, nearest_peptides=1, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> list[str]:
        """Suggest double mutant sequences for a given sequence to reduce PMB (avg distance to nearest reference peptides)

        NOTE that this does not guarantee complete humanness - but it should
        reduce the "PMB" - average distance to nearest human peptide.

        ALSO NOTE this does NOT take into account the frequency of the hit peptides (how many sequences they appear in)
        but rather the number of unique nearest peptides in the reference db that contain the given mutation.

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param allowed_amino_acids: Allowed amino acids for mutations
        :return: List of double mutant sequences, from most frequent
        """
        pssm_mutations = self.compute_pssm(seq, nearest_peptides=nearest_peptides, ignore_wildtype=True, allowed_amino_acids=allowed_amino_acids)
        # Flatten the DataFrame
        flat_cells = pssm_mutations.unstack()  # MultiIndex (column, row)
        # Filter out zeros
        flat_cells = flat_cells[flat_cells != 0]
        # Generate all combinations of cell pairs
        cell_pairs = itertools.combinations(flat_cells.items(), 2)
        # Compute products and sort, discard mutations at same position
        sorted_pairs = sorted(
            ((a_key, b_key, a_val * b_val) for (a_key, a_val), (b_key, b_val) in cell_pairs if a_key[1] != b_key[1]),
            key=lambda x: x[-1],
            reverse=True
        )
        mutants = []
        for (mut_aa1, pos1), (mut_aa2, pos2), product in sorted_pairs:
            mutant = list(seq)
            mutant[pos1 - 1] = mut_aa1
            mutant[pos2 - 1] = mut_aa2
            mutants.append(''.join(mutant))
        return mutants

    def find_double_mutant_seqs(self, seq: str, nearest_peptides=1, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY'):
        """Exhaustively find all possible double mutant sequences for a given sequence where all overlapping peptides
        are found in the reference database (are human).

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param allowed_amino_acids: Allowed amino acids for mutations
        :return: List of double mutant sequences (if any) whose peptides are contained in the db
        """
        candidates = self.suggest_double_mutant_candidates(seq, nearest_peptides=nearest_peptides, allowed_amino_acids=allowed_amino_acids)
        return [candidate for candidate in candidates if self.compute_peptide_content(candidate) == 1.0]

    def suggest_triple_mutant_candidates(self, seq: str, nearest_peptides=1, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY') -> list[str]:
        """Suggest triple mutant sequences for a given sequence to reduce PMB (avg distance to nearest reference peptides)

        NOTE that this does not guarantee complete humanness - but it should
        reduce the "PMB" - average distance to nearest human peptide.

        ALSO NOTE this does NOT take into account the frequency of the hit peptides (how many sequences they appear in)
        but rather the number of unique nearest peptides in the reference db that contain the given mutation.

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param allowed_amino_acids: Allowed amino acids for mutations
        :return: List of triple mutant sequences, from most frequent
        """
        pssm_mutations = self.compute_pssm(seq, nearest_peptides=nearest_peptides, ignore_wildtype=True, allowed_amino_acids=allowed_amino_acids)
        # Flatten the DataFrame
        flat_cells = pssm_mutations.unstack()  # MultiIndex (column, row)
        # Filter out zeros
        flat_cells = flat_cells[flat_cells != 0]
        # Generate all combinations of cell pairs
        cell_triples = itertools.combinations(flat_cells.items(), 3)
        # Compute products and sort, discard mutations at same position
        sorted_triples = sorted(
            ((a_key, b_key, c_key, a_val * b_val) for (a_key, a_val), (b_key, b_val), (c_key, c_val) in cell_triples if a_key[1] != b_key[1] and a_key[1] != c_key[1] and b_key[1] != c_key[1]),
            key=lambda x: x[-1],
            reverse=True
        )
        mutants = []
        for (mut_aa1, pos1), (mut_aa2, pos2), (mut_aa3, pos3), product in sorted_triples:
            mutant = list(seq)
            mutant[pos1 - 1] = mut_aa1
            mutant[pos2 - 1] = mut_aa2
            mutant[pos3 - 1] = mut_aa3
            mutants.append(''.join(mutant))
        return mutants

    def find_triple_mutant_seqs(self, seq: str, nearest_peptides=1, allowed_amino_acids='ACDEFGHIKLMNPQRSTVWY'):
        """Exhaustively find all possible triple mutant sequences for a given sequence where all overlapping peptides
        are found in the reference database (are human).

        :param seq: Sequence to analyze
        :param nearest_peptides: Number of nearest peptides to find, or list of lists of nearest peptides from db.find_nearest_peptides()
        :param allowed_amino_acids: Allowed amino acids for mutations
        :return: List of triple mutant sequences (if any) whose peptides are contained in the db
        """
        candidates = self.suggest_triple_mutant_candidates(seq, nearest_peptides=nearest_peptides, allowed_amino_acids=allowed_amino_acids)
        return [candidate for candidate in candidates if self.compute_peptide_content(candidate) == 1.0]


def add_wildcard(peptide: str, pos: int):
    aas = list(peptide)
    aas[pos] = '*'
    return ''.join(aas)

