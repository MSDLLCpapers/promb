from __future__ import annotations
import pandas as pd
import gzip
import sys
from typing import Iterable
import os
import glob
import re


RESIDUE_NAMES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 
}


def parse_pdb_atoms(pdb_path):
    df = pd.read_fwf(
        pdb_path,
        names=['record', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
            'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge'], 
        colspecs=[(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
            (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
            (78, 80)]
    )
    df = df[df['record'] == 'ATOM']
    return df


def iterate_sequences(input_paths: list[str], chains: str = None) -> Iterable[tuple[str, str, str]]:
    fasta_extensions = ('.fa', '.fasta', '.faa')
    text_extensions = ('.txt', '.txt.gz')
    pdb_extensions = ('.pdb', )
    extensions = fasta_extensions + text_extensions + pdb_extensions
    for i, input_path in enumerate(input_paths):
        filename = os.path.basename(input_path)
        print(f"Processing file {i+1}/{len(input_paths)}: {filename}", file=sys.stderr)
        if os.path.isdir(input_path):
            paths = [p for p in glob.glob(os.path.join(input_path, "*.*")) if p.endswith(extensions)]
            if not paths and len(input_paths) == 1:
                raise FileNotFoundError(f"No supported files found in directory: {input_path}")
            for result in iterate_sequences(paths, chains=chains):
                yield result
        elif input_path.lower().endswith(fasta_extensions):
            with open(input_path, 'r') as file:
                current_id = None
                current_sequence = []
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id is not None:
                            yield filename, current_id, ''.join(current_sequence)
                        current_id = line[1:].lstrip().split()[0]
                        current_sequence = []
                    else:
                        current_sequence.append(line)
                # Yield the last sequence
                if current_id is not None:
                    yield filename, current_id, ''.join(current_sequence)
        elif input_path.endswith(text_extensions):
            if input_path.endswith('.txt'):
                with open(input_path, 'r') as file:
                    yield filename, 'seq1', [line.strip().split()[0] for line in file.readlines()]
            elif input_path.endswith('.txt.gz'):
                with gzip.open(input_path, 'rt') as file:
                    yield filename, 'seq1', [line.strip().split()[0] for line in file.readlines()]
            else:
                raise NotImplementedError(f'Unsupported format: {input_path}')
        elif input_path.endswith(pdb_extensions):
            assert chains is not None, f'Please provide chains to parse when processing PDB files: {input_path}'
            assert isinstance(chains, str), f'Expected chains as single string, got: {chains}'
            atoms = parse_pdb_atoms(input_path)
            chains_list = list(chains.replace(',', ''))
            for chain in chains_list:
                residues = atoms[atoms['chainid'] == chain][['resseq','icode','resname']].drop_duplicates()
                seq = ''.join(residues['resname'].apply(lambda resname: RESIDUE_NAMES.get(resname, 'X')))
                yield filename, chain, seq
        else:
            raise NotImplementedError(f'Unsupported input file format: {input_path}, expected directory or one of {extensions}')



def parse_fasta(fasta_lines: Iterable[str]) -> dict[str, str]:
    """Parse fasta into a dictionary (header -> sequence)"""
    sequences = {}
    seq_id = None
    seq_lines = []

    assert not isinstance(fasta_lines, str), 'Expected list of lines or file handle, got string'
    
    for line in fasta_lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id:
                sequences[seq_id] = ''.join(seq_lines)
            seq_id = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    
    if seq_id:
        sequences[seq_id] = ''.join(seq_lines)

    for seq_id, seq in sequences.items():
        assert re.match(r'^[A-Z]+$', seq), f'Invalid sequence {seq_id}: {seq}'
    
    return sequences
