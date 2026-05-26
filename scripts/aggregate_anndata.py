#!/usr/bin/env python3
import argparse
import os
import re
import numpy as np
import pandas as pd
import anndata as ad

ROWS = list('ABCDEFGH')


def sample_num_to_well(n):
    '''S1=A1, S2=B1, ..., S8=H1, S9=A2, ..., S96=H12, S97 wraps to plate 2.'''
    n_within_plate = ((n - 1) % 96) + 1
    row_idx = (n_within_plate - 1) % 8
    col = (n_within_plate - 1) // 8 + 1
    return ROWS[row_idx] + str(col), n_within_plate


def extract_well_info(sample_name):
    match = re.search(r'[Ss](\d+)$', sample_name)
    if not match:
        print(f"[WARN] No S<n> token found in: {sample_name}")
        return {'BAM_name': sample_name, 'well_ID': None, 'well_nr': None,
                'sample_ID': None, 'seq_plate_ID': None}

    n = int(match.group(1))
    well_id, well_nr = sample_num_to_well(n)
    return {
        'BAM_name': sample_name,
        'well_ID': well_id,
        'well_nr': well_nr,
        'sample_ID': f'S{n}',
        'seq_plate_ID': str((n - 1) // 96 + 1),
    }


def main():
    parser = argparse.ArgumentParser(description='Aggregate per-sample recoding CSVs into an AnnData')
    parser.add_argument('--csvs', nargs='+', required=True, help='Per-sample CSVs (one per BAM)')
    parser.add_argument('--output', required=True, help='Output .h5ad path')
    args = parser.parse_args()

    csvs = sorted(args.csvs)
    print(f"Aggregating {len(csvs)} CSVs")

    samples = []
    dfs = []
    for csv in csvs:
        sample_name = os.path.basename(csv).replace('_recoding_analysis.csv', '')
        df = pd.read_csv(csv).sort_values('position').reset_index(drop=True)
        samples.append(sample_name)
        dfs.append(df)

    positions = dfs[0]['position'].to_list()
    for df, name in zip(dfs, samples):
        if df['position'].to_list() != positions:
            raise ValueError(f"Position vector mismatch in {name}; all CSVs must share the same positions")

    X_frequency = np.array([df['frequency'].to_numpy() for df in dfs], dtype=float)
    X_depth     = np.array([df['depth'].to_numpy()     for df in dfs], dtype=int)
    X_recoded   = np.array([df['recoded_codons'].to_numpy() for df in dfs], dtype=int)

    obs = pd.DataFrame([extract_well_info(s) for s in samples])
    obs.index = samples

    var = pd.DataFrame({'position': positions})
    var.index = [f'codon_{p}' for p in positions]

    adata = ad.AnnData(
        X=X_frequency,
        obs=obs,
        var=var,
        layers={'depth': X_depth, 'recoded_count': X_recoded},
    )

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}  shape={adata.shape}")


if __name__ == '__main__':
    main()
