from __future__ import annotations
import argparse
import sys
import pandas as pd

from promb import init_db
from tqdm import tqdm

from promb.netmhc import predict_netmhc_top_ranks
from promb.parsing import iterate_sequences

def run_promb(
        db,
        input_paths,
        netmhc_alleles: list[str] = None,
        netmhc_peptide_lengths=None,
        mode='mutations',
        chains: str = None,
    ):
    if netmhc_peptide_lengths is None:
        netmhc_peptide_lengths = [12, 13, 14, 15, 16, 17]
    tuples = iterate_sequences(input_paths, chains=chains)
    skipped_ids = []
    dfs = []
    for filename, sequence_id, sequence in tqdm(tuples):
        if len(sequence) < db.peptide_length:
            skipped_ids.append(sequence_id)
            continue
        peptides = db.chop_seq_peptides(sequence)
        df = pd.DataFrame([{
            'File': filename,
            'ID': sequence_id,
            'Position': position,
            'Peptide': peptide,
        } for position, peptide in enumerate(peptides, start=1)])
        if mode == 'nearest':
            nearest = db.find_nearest_peptides(peptides, n=1)
            df['Nearest'] = [n[0] for n in nearest]
            df['Mutations'] = [sum(aa != bb for aa, bb in zip(peptide, n[0])) for peptide, n in zip(peptides, nearest)]
        elif mode == 'exact':
            df['Found'] = [db.contains(peptide) for peptide in peptides]
        elif mode == 'content':
            df['Content'] = [db.contains(peptide) for peptide in peptides]
            df = df.groupby(['File', 'ID'], as_index=False).agg({'Content': 'mean'})
            print(df)
        else:
            raise ValueError(f'Unknown mode: {mode}')
        if netmhc_alleles:
            core_ranks, core_binders, cores = predict_netmhc_top_ranks(
                sequence,
                k=db.peptide_length,
                alleles=netmhc_alleles,
                peptide_lengths=netmhc_peptide_lengths
            )
            df['Top_Rank'] = core_ranks.min(axis=1).values
            df['Top_Rank_Allele'] = core_ranks.idxmin(axis=1).values
            df['Top_Rank_Core'] = [cores.iloc[i][allele] for i, allele in enumerate(df['Top_Rank_Allele'])]
            assert (df['Peptide'] == core_ranks.index).all(), 'Expected the same order of peptides in df and core_ranks'
            df = pd.concat([
                df.reset_index(drop=True),
                core_ranks.add_suffix('_Rank').reset_index(drop=True),
                cores.add_suffix('_Core').reset_index(drop=True),
                core_binders.add_suffix('_Binder').reset_index(drop=True)
            ], axis=1)
        dfs.append(df)
    if not dfs:
        raise ValueError(f'No input sequences found in: {input_paths}')
    result = pd.concat(dfs)
    return result, skipped_ids


def run_and_save(db_path: str, peptide_length=None, output_path=None, **kwargs):
    db = init_db(db_path, peptide_length=peptide_length)
    result, skipped_ids = run_promb(db=db, **kwargs)

    if output_path and output_path != '-':
        if output_path.endswith(('.csv', '.csv.gz')):
            result.to_csv(output_path, index=False)
        elif output_path.endswith('.xlsx'):
            with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
                result.to_excel(writer, sheet_name='Peptides', index=False)
                sheet = writer.sheets['Peptides']
                header = result.columns.tolist()
                id_col = header.index('ID')
                sheet.set_column(id_col, id_col, 20)
                peptide_col = header.index('Peptide')
                sheet.set_column(peptide_col, peptide_col, db.peptide_length + 2)
                # format Mutations from white (0) to red (3)
                if 'Mutations' in header:
                    mutations_col = header.index('Mutations')
                    sheet.conditional_format(1, mutations_col, len(result), mutations_col, {
                        'type': '2_color_scale',
                        'min_type': 'num',
                        'min_value': 0,
                        'max_type': 'num',
                        'max_value': 4,
                        'min_color': '#FFFFFF',
                        'max_color': '#ff7b91'
                    })
                # format Rank columns from Dark purple (1.0) to purple (5.0) to white (20.0)
                if 'Top_Rank' in header:
                    rank_from = header.index('Top_Rank')
                    rank_to = [(i, col) for i, col in enumerate(header) if col.endswith('_Rank')][-1][0]
                    sheet.conditional_format(1, rank_from, len(result), rank_to, {
                        'type': '3_color_scale',
                        'min_type': 'num',
                        'min_value': 1,
                        'mid_type': 'num',
                        'mid_value': 5,
                        'max_type': 'num',
                        'max_value': 20,
                        'min_color': '#c362ff',
                        'mid_color': '#cb98ff',
                        'max_color': '#FFFFFF'
                    })
        elif output_path == '/dev/null':
            pass
        else:
            raise NotImplementedError(f'Unsupported output format, expected .csv or .xlsx, got: {output_path}')
    else:
        result.to_csv(sys.stdout, index=False)

    numeric_columns = [c for c in result.select_dtypes(include=['number', 'bool']).columns if c not in ['Position']]
    for column in numeric_columns:
        print(f"{column}: Statistics across sequences", file=sys.stderr)
        print(result.groupby(['File', 'ID'])[column].mean().describe().T, file=sys.stderr)

    if skipped_ids:
        print(f"NOTE: Skipped {len(skipped_ids)} sequences because they were shorter "
              f"than the peptide length ({db.peptide_length})", file=sys.stderr)
    if output_path:
        print("Saved to:", output_path, file=sys.stderr)


def run_oasis(input_paths, output_path=None, chains=None):
    db = init_db('human-oas')
    tuples = iterate_sequences(input_paths, chains=chains)
    
    result = []
    skipped_ids = []
    for filename, sequence_id, sequence in tqdm(tuples):
        result.append({
            'File': filename,
            'ID': sequence_id,
            'Identity': db.compute_peptide_content(sequence),
        })
    result = pd.DataFrame(result)

    if output_path and output_path != '-':
        if output_path.endswith(('.csv', '.csv.gz')):
            result.to_csv(output_path, index=False)
        elif output_path.endswith('.xlsx'):
            with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
                result.to_excel(writer, sheet_name='Scores', index=False)
                sheet = writer.sheets['Scores']
                header = result.columns.tolist()
                id_col = header.index('ID')
                sheet.set_column(id_col, id_col, 20)
                identity_col = header.index('Identity')
                sheet.conditional_format(1, identity_col, len(result), identity_col, {
                    'type': '2_color_scale',
                    'min_type': 'num',
                    'min_value': 0,
                    'max_type': 'num',
                    'max_value': 1,
                    'min_color': '#ff7b91',
                    'max_color': '#FFFFFF'
                })
        elif output_path == '/dev/null':
            pass
        else:
            raise NotImplementedError(f'Unsupported output format, expected .csv or .xlsx, got: {output_path}')
    else:
        result.to_csv(sys.stdout, index=False)

    if skipped_ids:
        print(f"NOTE: Skipped {len(skipped_ids)} sequences because they were shorter "
              f"than the peptide length ({db.peptide_length})", file=sys.stderr)
    if output_path:
        print("Saved to:", output_path, file=sys.stderr)

def main(args=None):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')

    nearest_parser = subparsers.add_parser('nearest')
    exact_parser = subparsers.add_parser('exact')
    oasis_parser = subparsers.add_parser('oasis')
    
    for p in [nearest_parser, exact_parser, oasis_parser]:
        p.add_argument('input', nargs='+', help='Input file(s) - FASTA, .txt/.txt.gz (with one peptide per line), .pdb file, or a directory of such files')
        p.add_argument('-c', '--chains', help='Chains to parse when processing PDB files (comma-separated)')
        p.add_argument('-o', '--output', help='Output CSV file, will print to stdout if not provided')

    for p in [nearest_parser, exact_parser]:
        p.add_argument('-d', '--db', required=True, help='Input database name (oas-human, oas-swissprot) or path to .txt.gz file (one peptide per line, same length) or .fasta.gz file (whole sequences in fasta format)')
        p.add_argument('-l', '--peptide-length', type=int, help='Length of peptide (when using oas-swissprot or other fasta database)')
        p.add_argument('-a', '--netmhc-alleles', help='Run netMHCIIpan for given alleles (comma-separated, for example DRB1_0101,DRB1_0102). '
                                                            'requires netMHCIIpan license and netMHCIIpan to be installed and on PATH')

    args = parser.parse_args(args)

    if args.command is None:
        parser.print_help()
    if args.command == 'oasis':
        run_oasis(
            input_paths=args.input,
            chains=args.chains,
            output_path=args.output,
        )
    else:
        run_and_save(
            db_path=args.db,
            input_paths=args.input,
            peptide_length=args.peptide_length,
            chains=args.chains,
            netmhc_alleles=args.netmhc_alleles.split(',') if args.netmhc_alleles else None,
            mode=args.command,
            output_path=args.output,
        )


if __name__ == '__main__':
    main()
