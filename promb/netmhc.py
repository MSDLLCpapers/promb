from __future__ import annotations
import os
import subprocess
import sys
import tempfile
from io import StringIO

import numpy as np
import pandas as pd


# NetMHCIIpan main arguments
#
# 	PARAMETER            DEFAULT VALUE        DESCRIPTION
# 	[-f filename]                             File with the input data
# 	[-BA]                0                    Include BA predictions, default is EL only
# 	[-context]           0                    Predict with context encoding
# 	[-hlapseudo filename] $NETMHCIIpan/data/pseudosequence.2023.all.X.dat File with HLA pseudo sequences
# 	[-a filename]        DRB1_0101            Allele for prediction (multiple alleles with ,)
# 	[-choose]            0                    Choose alpha and beta chains separately
# 	[-cha string]                             Alpha chain name
# 	[-chb string]                             Beta chain name
# 	[-rankS float]       1.000000             Threshold for strong binders (%Rank)
# 	[-rankW float]       5.000000             Threshold for weak binders (%Rank)
# 	[-filter]            0                    Toggle filtering of output
# 	[-rankF float]       10.000000            Threshold for filtering output (%Rank), if -filter option in on
# 	[-hlaseq filename]                        File with full length MHC beta chain sequence (used alone for HLA-DR)
# 	[-hlaseqA filename]                       File with full length MHC alpha chain sequence (used with -hlaseq option)
# 	[-inptype int]       0                    Input type [0] FASTA [1] Peptide
# 	[-length filename]   15                   Peptide length (multiple length with ,). Used for FASTA input only.
# 	[-s]                 0                    Sort output on descending affinity
# 	[-u]                 0                    Print unique binding core only
# 	[-xls]               0                    Save output into xls file
# 	[-xlsfile filename]  NetMHCIIpan_out.xls  Filename for xls dump
# 	[-list]              0                    Print the list of possible alleles and exit
# 	[-termAcon]          0                    Encode context extending protein sequence as A [default is X]
# 	[-inv_all]           0                    Consider inversion for all molecules
# 	[-v]                 0                    Verbose mode


def predict_netmhc_top_ranks(sequence: str, k: int, alleles: list[str], peptide_lengths: list[int]):
    """For each overlapping k-mer (of fixed length) in a sequence, predict the top ranking MHC binder peptides containing that k-mer"""
    peptides = [sequence[i: i + k] for i in range(len(sequence) - k + 1)]
    df = run_netmhc(
        sequence=sequence,
        alleles=alleles,
        peptide_lengths=peptide_lengths
    )
    binders = []
    for allele in alleles:
        allele_df = df[allele].copy()
        allele_df.insert(0, 'Allele', allele)
        allele_df.insert(1, 'StartPos', df[('Input', 'Pos')])
        allele_df.insert(2, 'EndPos', df[('Input', 'Pos')] + df[('Input', 'Peptide')].str.len() - 1)
        allele_df.insert(3, 'Binder', df[('Input', 'Peptide')])
        binders.append(allele_df)
    binders = pd.concat(binders).sort_values(by=['Allele', 'StartPos', 'EndPos'])
    # Allele     StartPos  EndPos  Binder             Core       Inverted  Score     Rank
    # DRB1_0101  1         12      EIVLTQSPATLS       VLTQSPATL  0         0.028873  16.26001
    # DRB1_0101  1         13      EIVLTQSPATLSL      VLTQSPATL  0         0.024097  17.480165
    # DRB1_0101  1         14      EIVLTQSPATLSLS     VLTQSPATL  0         0.044627  13.456401
    # DRB1_0101  1         15      EIVLTQSPATLSLSP    LTQSPATLS  0         0.059213  11.800208
    # DRB1_0101  1         16      EIVLTQSPATLSLSPG   LTQSPATLS  0         0.027442  16.597033
    # DRB1_0101  1         17      EIVLTQSPATLSLSPGE  LTQSPATLS  0         0.012217  22.637009
    # DRB1_0101  2         13      IVLTQSPATLSL       VLTQSPATL  0         0.002928  35.898647
    # DRB1_0101  2         14      IVLTQSPATLSLS      VLTQSPATL  0         0.003531  33.989304
    # DRB1_0101  2         15      IVLTQSPATLSLSP     LTQSPATLS  0         0.0072    27.152369
    # DRB1_0101  2         16      IVLTQSPATLSLSPG    LTQSPATLS  0         0.005733  29.296297
    # DRB1_0101  2         17      IVLTQSPATLSLSPGE   LTQSPATLS  0         0.001832  41.084415
    # DRB1_0101  2         18      IVLTQSPATLSLSPGER  LTQSPATLS  0         0.000753  51.862068
    top_ranks = []
    top_binders = []
    top_cores = []
    for pos, peptide in enumerate(peptides, start=1):
        binders_with_peptide = binders.query('StartPos <= @pos and EndPos >= (@pos + @k - 1)')
        assert binders_with_peptide['Binder'].str.contains(peptide).all(), f'Peptide {peptide} not present in expected positions'
        top_by_allele = binders_with_peptide.groupby('Allele', sort=True).apply(lambda rows: rows.sort_values(by='Rank').iloc[0])
        # Allele     StartPos  EndPos  Binder             Core       Inverted  Score     Rank
        # DRB1_0101  1         15      EIVLTQSPATLSLSP    LTQSPATLS  0         0.059213  11.800208
        # DRB1_0102  1         14      EIVLTQSPATLSLS     TQSPATLSL  0         0.099213  10.800208
        top_ranks.append(top_by_allele['Rank'].rename(peptide))
        top_binders.append(top_by_allele['Binder'].rename(peptide))
        top_cores.append(top_by_allele['Core'].rename(peptide))
    top_ranks = pd.DataFrame(top_ranks)
    top_ranks.columns.name = None
    top_ranks.index.name = 'Peptide'
    # Peptide	DRB1_0101   DRB1_0102
    # CQHSRDLPL 98.5        99.0
    # CRASKGVST 58.1        60.0
    top_binders = pd.DataFrame(top_binders)
    top_binders.columns.name = None
    top_binders.index.name = 'Peptide'
    # Peptide	DRB1_0101   DRB1_0102
    # CQHSRDLPL YCQHSRDLPLTFG YCQHSRDLPLTF
    # CRASKGVST TLSCRASKGVSTSGY ATLSCRASKGVSTSG
    top_cores = pd.DataFrame(top_cores)
    top_cores.columns.name = None
    top_cores.index.name = 'Peptide'
    # Peptide	DRB1_0101   DRB1_0102
    # CQHSRDLPL YCQHSRDLP   HSRDLPLTF
    # CRASKGVST TLSCRASKG   TLSCRASKG
    assert top_ranks.index.equals(top_binders.index), 'Expected the same peptides in top_ranks and top_binders'
    return top_ranks, top_binders, top_cores


def run_netmhc(sequence: str, alleles: list[str], peptide_lengths: list[int]) -> pd.DataFrame:
    """
    Run NetMHC for a protein sequence.

    :param sequence: Protein sequence
    :param alleles: List of alleles
    :param peptide_lengths: List of peptide lengths
    :return: DataFrame with NetMHC output
    """
    # sort alleles because netMHC will also sort them in the output table
    alleles = sorted(alleles)
    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, 'input.fa')
        output_path = os.path.join(temp_dir, 'output.tsv')
        with open(input_path, 'wt') as f:
            f.write('>sequence\n')
            f.write(sequence + '\n')
        r = subprocess.run([
            'netMHCIIpan',
            '-f', input_path,
            '-inptype', '0',
            '-xls', '-xlsfile', output_path,
            '-length', ','.join(map(str, peptide_lengths)),
            '-a', ','.join(alleles)
        ], capture_output=True, text=True)
        if r.returncode or not os.path.exists(output_path):
            print(r.stdout, file=sys.stderr)
            print(r.stderr, file=sys.stderr)
            raise Exception("NetMHCIIpan failed, see error above")
        io = StringIO()
        with open(output_path) as f:
            output_alleles = [allele for allele in f.readline().strip().split('\t') if allele]
            assert output_alleles == alleles, 'Expected alleles {} in output, got {}'.format(alleles, output_alleles)
            second_line = f.readline().strip().split('\t')
            assert second_line == ['Pos', 'Peptide', 'ID', 'Target'] + ['Core', 'Inverted', 'Score', 'Rank'] * len(alleles) + ['Ave', 'NB'], \
                'Unexpected NetMHCIIpan output format: ' + str(second_line)
            fixed_header = ['Input', 'Input', 'Input', 'Input'] + [allele for allele in alleles for _ in range(4)] + ['Summary', 'Summary']
            io.write('\t'.join(fixed_header))
            io.write('\n')
            io.write('\t'.join(second_line))
            io.write('\n')
            io.write(f.read())
        io.seek(0)
        return pd.read_csv(io, sep='\t', header=[0, 1])
