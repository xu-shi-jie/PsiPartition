import os
import sys
from functools import partial
from pathlib import Path

sys.path.append('.')
import argparse

import numpy as np
import pandas as pd
from loguru import logger

import wandb
from PsiPartition_wandb import calculate_sorting
from utils import (default_iqtree, dna_str, extract_stat, remove_caches,
                   write_part_file)


def PsiPartition(msa: str, format: str, alphabet: str, k: int, w: dict = None):
    _, sorting = calculate_sorting(msa, format, alphabet, w)
    indices = np.digitize(sorting, np.linspace(
        sorting.min(), sorting.max(), k))
    part_file = Path(msa).with_suffix('.parts')
    write_part_file(part_file, indices)

    with open(Path(msa).with_suffix('.log'), 'w') as log_file:
        remove_caches(Path(msa))
        log_file = open(Path(msa).with_suffix('.log'), 'w')
        default_iqtree([
            '-s', msa, '-pre', Path(msa).with_suffix(''),
            '-spp', Path(msa).with_suffix('.parts'), '-nt', '8',],
            log_file)
    stat = extract_stat(Path(msa).with_suffix('.iqtree'))
    return stat

def opt_func(args) -> float:
    wandb.init(project='PsiPartition')
    msa = args.msa
    format = args.format
    alphabet = args.alphabet
    k = wandb.config.k

    w = {
        'A': wandb.config.A, 'C': wandb.config.C,
        'G': wandb.config.G, 'T': wandb.config.T, '-': wandb.config.UNK, }
    s = sum(w.values())
    w = {k: v/s for k, v in w.items()}

    stat = PsiPartition(msa, format, alphabet, k, w)
    wandb.log(stat)
    # write results to file
    with open(Path(msa).with_suffix(f'.{args.opt}.csv'), 'a') as f:
        f.write(
            f'{k},{stat["BIC"]},{stat["AICc"]},{stat["log-likelihood"]},{stat["num_params"]},')
        if w is not None:
            f.write(
                ','.join([f'{w[k]}' for k in dna_str]))
        f.write('\n')

    # move iqtree, treefile, part file to PsiPartition directory
    sweep_id = wandb.run.name.split('-')[-1]
    iter_dir = Path(msa).parent / str(Path(msa).stem) / str(args.opt) / f'Iter{sweep_id}'
    logger.info(f'Moving results to {iter_dir}')
    iter_dir.mkdir(exist_ok=True, parents=True)
    os.system(f'mv {Path(msa).with_suffix(".iqtree")} {iter_dir}')
    os.system(f'mv {Path(msa).with_suffix(".treefile")} {iter_dir}')
    os.system(f'mv {Path(msa).with_suffix(".parts")} {iter_dir}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa', type=str, default='data/BayesOpt/Morpho.fasta')
    parser.add_argument('--format', type=str, default='fasta')
    parser.add_argument('--alphabet', type=str, default='dna')
    parser.add_argument('--n_iter', type=int, default=60)
    parser.add_argument('--opt', type=str, default='bayes')
    args = parser.parse_args()
    
    n_iter = args.n_iter
    file = Path(args.msa)
    ########################## test DNA empirical ##########################
    part_file = file.with_suffix('.parts')
    csv_file = Path(file).with_suffix(f'.{args.opt}.csv')
    if csv_file.exists():
        df = pd.read_csv(csv_file, header=None)
        if df.shape[0] != n_iter:
            csv_file.unlink()
            logger.info(f'Records {df.shape[0]} less than {n_iter}, will rerun PsiPartition')
        else:
            logger.info(f'Found existing PsiPartition results for {file}')
    else:
        logger.info(f'CSV file {csv_file} does not exist, running PsiPartition')

    if not csv_file.exists():
        # run PsiPartition
        logger.info(f'Running PsiPartition on {file}')
        params = {
            'k': {'min': 1, 'max': 30, 'distribution': 'int_uniform'},
            'A': {'min': 0, 'max': 1, 'distribution': 'uniform'},
            'C': {'min': 0, 'max': 1, 'distribution': 'uniform'},
            'G': {'min': 0, 'max': 1, 'distribution': 'uniform'},
            'T': {'min': 0, 'max': 1, 'distribution': 'uniform'},
            'UNK': {'min': 0, 'max': 1, 'distribution': 'uniform'},
        }
        sweep_configuration = {
            'method': args.opt,
            "name": "PsiPartition",
            'metric': {'name': 'BIC', 'goal': 'minimize'}, 
            'parameters': params,}
        sweep_id = wandb.sweep(sweep_configuration, project='PsiPartition')
        wandb.agent(sweep_id, partial(opt_func, args), count=n_iter)

    # check csv file and obtain the best parameters
    df = pd.read_csv(csv_file, header=None)
    best = df.iloc[df[1].idxmin()]
    k = int(best[0])
    w = {
        'A': best[5], 'C': best[6], 'G': best[7], 'T': best[8], '-': best[9], }
    s = sum(w.values())
    w = {k: v/s for k, v in w.items()}
    # run IQ-TREE
    logger.info(f'Running IQ-TREE on {file} with k={k}, w={w}')
    PsiPartition(file, 'fasta', 'dna', k, w)
