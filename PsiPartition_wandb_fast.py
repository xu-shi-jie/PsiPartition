"""
A script to evaluate the sorting partitioning method.
@Author: Shijie Xu
@Date: 2024-01-02
"""
import argparse
from collections import defaultdict
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import wandb
from Bio import AlignIO
from loguru import logger
from tqdm import tqdm

from utils import (aa_str, default_iqtree, dna_str, extract_stat,
                   remove_caches, write_part_file)


def calculate_sorting(file_path: str, format, datatype, w) -> np.ndarray:
    msa = [str(r.seq) for r in AlignIO.read(file_path, format)]
    n_seqs, n_sites = len(msa), len(msa[0])
    alphabet = aa_str if datatype == 'aa' else dna_str

    site_partitions = []
    for j, s in enumerate(zip(*msa)):
        partition = defaultdict(list)
        for i, c in enumerate(s):
            if c in alphabet[:-1]:
                partition[c].append(i)
            else:
                partition['-'].append(i)
        site_partitions.append({
            k: frozenset(v) for k, v in partition.items()})
    corr = np.ones((n_sites, n_sites))
    for i in tqdm(range(n_sites), desc="Computing sorting matrix", leave=False):
        for j in range(n_sites):
            if i != j:
                s = 0
                for c2, p2 in site_partitions[j].items():
                    for c1, p1 in site_partitions[i].items():
                        if p2.issubset(p1):
                            s += w[c2] if w else 1
                            break
                corr[i, j] = s / len(site_partitions[j]) if s else 0
    sorting = (corr < corr.T).sum(1)

    return corr, sorting


def PsiPartition(msa: str, format: str, alphabet: str, k: int, w: dict = None):
    _, sorting = calculate_sorting(msa, format, alphabet, w)
    log_file = open(Path(msa).with_suffix('.log'), 'w')
    indices = np.digitize(sorting, np.linspace(
        sorting.min(), sorting.max(), k))
    part_file = Path(msa).with_suffix('.parts')
    write_part_file(part_file, indices)
    with open(Path(msa).with_suffix('.log'), 'w') as log_file:
        remove_caches(Path(msa))
        default_iqtree([
            '-s', msa, '-pre', Path(msa).with_suffix(''),
            '-spp', Path(msa).with_suffix('.parts'), '-nt', '4',],
            log_file)
    stat = extract_stat(Path(msa).with_suffix('.iqtree'))
    return stat


def opt_func(args) -> float:
    wandb.init(project='PsiPartition')
    msa = args.msa
    format = args.format
    alphabet = args.alphabet
    k = wandb.config.k
    w = None

    stat = PsiPartition(msa, format, alphabet, k, w)

    wandb.log(stat)
    # write results to file
    with open(Path(msa).with_suffix('.csv'), 'a') as f:
        f.write(
            f'{k},{stat["BIC"]},{stat["AICc"]},{stat["log-likelihood"]},{stat["num_params"]},')
        if w is not None:
            f.write(
                ','.join([f'{w[k]}' for k in (aa_str if alphabet == "aa" else dna_str)]))
        f.write('\n')

if __name__ == '__main__':
    # fmt: off
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa', type=str, default='data/DNA/empirical/Morpho.fasta', help='MSA file')
    parser.add_argument('--format', type=str, default='fasta', help='fasta or phylip')
    parser.add_argument('--alphabet', type=str, default='dna', help='dna or aa')
    parser.add_argument('--max_partitions', type=int, default=30)
    parser.add_argument('--n_iter', type=int, default=10)
    args = parser.parse_args()
    # fmt: on

    tuning_csv = Path(args.msa).with_suffix('.csv')
    if not tuning_csv.exists():
        # define bayesian optimization search space
        sweep_configuration = {
            "method": "bayes",
            "name": "PsiPartition",
            "metric": {"goal": "minimize", "name": "BIC"},
            "parameters": {
                "k": {"values": list(range(2, args.max_partitions+1))},
            },
        }
        sweep_id = wandb.sweep(sweep_configuration, project='SitePartition')
        wandb.agent(sweep_id, partial(opt_func, args), count=args.n_iter)

    # read tuning results
    df = pd.read_csv(tuning_csv, header=None)
    best = df.iloc[df[1].idxmin()]
    k = int(best[0])
    # run PsiPartition
    logger.info(f'Running PsiPartition on {args.msa}')
    PsiPartition(args.msa, args.format, args.alphabet, k, None)