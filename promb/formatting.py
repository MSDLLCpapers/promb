from __future__ import annotations
import sys

def chunk_list(lst: list, n: int) -> list[list]:
    """Return list of successive n-sized batches from lst."""
    return [lst[i : i + n] for i in range(0, len(lst), n)]


def print_nearest(queries: list[str], hit_lists: list[list[str]], wrap=10, file=sys.stdout):
    assert len(queries) == len(hit_lists), f"Expected same number of positions in query list and hit list, got {len(queries)} != {len(hit_lists)}"
    if not queries:
        return
    
    print("^ = mutation, + = similar residue", file=file)
    print(file=file)
    
    pos = 0
    for query_batch, hits_batch in zip(chunk_list(queries, wrap), chunk_list(hit_lists, wrap)):
        number_line = []
        query_line = []
        hit_lines = [[] for _ in range(len(hit_lists[0]))]
        mut_lines = [[] for _ in range(len(hit_lists[0]))]

        for query, hits in zip(query_batch, hits_batch):
            pos += 1
            number_line.append(str(pos).ljust(len(query)))
            query_line.append(query)
            for hit, hit_line, mut_line in zip(hits, hit_lines, mut_lines):
                hit_line.append(hit if hit != query else 'matches'.ljust(len(query)))
                mut_line.append(''.join(('+' if are_similar(aa, bb) else '^') if aa != bb else ' ' for aa, bb in zip(query, hit)))

        print('  '.join(number_line), file=file)
        print('  '.join(query_line) + ' queries', file=file)
        for i, (hit_line, mut_line) in enumerate(zip(hit_lines, mut_lines)):
            print('  '.join(hit_line) + (f' #{i+1} hit' if len(hit_lines) > 1 else ' hits'), file=file)
            print('  '.join(mut_line), file=file)
        print(file=file)

        
# These are the amino acid pairs that are considered similar when calculating sequence similarity.
# Based on positive scores in BLOSUM62 (should correspond to BLAST similarity)
SIMILAR_PAIRS = frozenset({
    'AA', 'AS', 'CC', 'DD', 'DE', 'DN', 'ED', 'EE', 'EK', 'EQ', 'FF', 'FW', 'FY', 'GG', 'HH', 'HN', 'HY',
    'II', 'IL', 'IM', 'IV', 'KE', 'KK', 'KQ', 'KR', 'LI', 'LL', 'LM', 'LV', 'MI', 'ML', 'MM', 'MV', 'ND',
    'NH', 'NN', 'NS', 'PP', 'QE', 'QK', 'QQ', 'QR', 'RK', 'RQ', 'RR', 'SA', 'SN', 'SS', 'ST', 'TS', 'TT',
    'VI', 'VL', 'VM', 'VV', 'WF', 'WW', 'WY', 'YF', 'YH', 'YW', 'YY',

    # ambiguous aminoacid codes
    'XX',
    'BB', 'BN', 'NB', 'BD', 'DB',
    'ZZ', 'ZQ', 'QZ', 'ZE', 'EZ',
    'JJ', 'JL', 'LJ', 'JI', 'IJ',
    'UU',
    'OO'
})

def are_similar(aa1: str, aa2: str) -> bool:
    """Check if amino acid pair is considered similar (have positive BLOSUM score)"""
    return aa1+aa2 in SIMILAR_PAIRS
