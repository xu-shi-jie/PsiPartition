import os
import sys
from pathlib import Path

sys.path.append('.')
import argparse

import pandas as pd
from loguru import logger

from PsiPartition_wandb import PsiPartition
from utils import default_iqtree, remove_caches

if __name__ == '__main__':
    psipartition_path = r'PsiPartition_wandb.py'

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=str, default='DNA-simulated-seqlen')
    parser.add_argument('--dir', type=str, default='')
    args = parser.parse_args()
    
    if args.dataset == 'DNA-empirical':
        ########################## test DNA empirical ##########################
        dna_empirical_dir = r'data/DNA/empirical'
        sp_dir = Path(dna_empirical_dir) / 'PsiPartition'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_empirical_dir).glob('*.fasta'):
            part_file = sp_dir / f'{file.stem}.parts'
            new_file = part_file.with_suffix(".fasta")
            # move file to PsiPartition directory
            os.system(f'cp {file} {new_file}')
            file = new_file
            
            csv_file = Path(file).with_suffix('.csv')
            if csv_file.exists():
                df = pd.read_csv(csv_file, header=None)
                if df.shape[0] != 30:
                    csv_file.unlink()
                    logger.info(f'Records {df.shape[0]} less than 30, will rerun PsiPartition')
                else:
                    logger.info(f'Found existing PsiPartition results for {file}')
            else:
                logger.info(f'CSV file {csv_file} does not exist, running PsiPartition')

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {file} --format fasta --alphabet dna --max_partitions 30 --n_iter 30')

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


    elif args.dataset == 'DNA-simulated':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated'
        sp_dir = Path(dna_simulated_dir) / 'PsiPartition'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_simulated_dir).glob(f'{args.dir}/*.fasta'):
            part_file = sp_dir / file.parent.name / f'{file.stem}.parts'
            part_file.parent.mkdir(exist_ok=True, parents=True)
            new_file =part_file.with_suffix(".fasta")
            os.system(f'cp {file} {new_file}')
            file = new_file

            csv_file = Path(file).with_suffix('.csv')
            if csv_file.exists():
                df = pd.read_csv(csv_file, header=None)
                if df.shape[0] != 30:
                    csv_file.unlink()
                    logger.info(f'Records {df.shape[0]} less than 30, will rerun PsiPartition')
                else:
                    logger.info(f'Found existing PsiPartition results for {file}')
            else:
                logger.info(f'CSV file {csv_file} does not exist, running PsiPartition')

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {new_file} --format fasta --alphabet dna --max_partitions 30 --n_iter 30')

            df = pd.read_csv(csv_file, header=None)
            best = df.iloc[df[1].idxmin()]
            k = int(best[0])
            w = {
                'A': best[5], 'C': best[6], 'G': best[7], 'T': best[8], '-': best[9], }
            s = sum(w.values())
            w = {k: v/s for k, v in w.items()}
            logger.info(f'Running IQ-TREE on {file} with k={k}, w={w}')
            PsiPartition(file, 'fasta', 'dna', k, w)


    elif args.dataset == 'protein-empirical':
        ########################## test protein empirical ##########################
        n_iter = 100
        n_part = 20
        protein_empirical_dir = r'data/Protein/empirical'
        sp_dir = Path(protein_empirical_dir) / 'PsiPartition'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(protein_empirical_dir).glob('Nguyen_2016e.fasta'):
            part_file = sp_dir / f'{file.stem}.parts'
            new_file = part_file.with_suffix(".fasta")
            os.system(f'cp {file} {new_file}')
            file = new_file

            csv_file = Path(file).with_suffix('.csv')
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
                os.system(f'python {psipartition_path} --msa {file} --format fasta --alphabet aa --max_partitions {n_part} --n_iter {n_iter}')

            df = pd.read_csv(csv_file, header=None)
            best = df.iloc[df[1].idxmin()]
            k = int(best[0])
            w = {
                'A': best[5], 'R': best[6], 'N': best[7], 'D': best[8], 'C': best[9],
                'Q': best[10], 'E': best[11], 'G': best[12], 'H': best[13], 'I': best[14],
                'L': best[15], 'K': best[16], 'M': best[17], 'F': best[18], 'P': best[19],
                'S': best[20], 'T': best[21], 'W': best[22], 'Y': best[23], 'V': best[24], 
                '-': best[25], }
            
            s = sum(w.values())
            w = {k: v/s for k, v in w.items()}
            logger.info(f'Running IQ-TREE on {file} with k = {k}, w = {w}')
            PsiPartition(file, 'fasta', 'aa', k, w)

    elif args.dataset == 'DNA-simulated-seqlen':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated-seqlen'
        sp_dir = Path(dna_simulated_dir) / 'PsiPartition'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_simulated_dir).glob(f'{args.dir}/*.fasta'):
            part_file = sp_dir / file.parent.name / f'{file.stem}.parts'
            part_file.parent.mkdir(exist_ok=True, parents=True)
            new_file =part_file.with_suffix(".fasta")
            os.system(f'cp {file} {new_file}')
            file = new_file

            csv_file = Path(file).with_suffix('.csv')
            if csv_file.exists():
                df = pd.read_csv(csv_file, header=None)
                if df.shape[0] != 30:
                    csv_file.unlink()
                    logger.info(f'Records {df.shape[0]} less than 30, will rerun PsiPartition')
                else:
                    logger.info(f'Found existing PsiPartition results for {file}')
            else:
                logger.info(f'CSV file {csv_file} does not exist, running PsiPartition')

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {new_file} --format fasta --alphabet dna --max_partitions 30 --n_iter 30')

            df = pd.read_csv(csv_file, header=None)
            best = df.iloc[df[1].idxmin()]
            k = int(best[0])
            w = {
                'A': best[5], 'C': best[6], 'G': best[7], 'T': best[8], '-': best[9], }
            s = sum(w.values())
            w = {k: v/s for k, v in w.items()}
            logger.info(f'Running IQ-TREE on {file} with k={k}, w={w}')
            PsiPartition(file, 'fasta', 'dna', k, w)