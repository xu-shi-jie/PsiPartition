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
    psipartition_path = r'PsiPartition_wandb_fast.py'

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=str, default='DNA-simulated-seqlen')
    args = parser.parse_args()
    
    if args.dataset == 'DNA-empirical':
        ########################## test DNA empirical ##########################
        dna_empirical_dir = r'data/DNA/empirical'
        sp_dir = Path(dna_empirical_dir) / 'PsiPartitionFast'
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
                if df.shape[0] != 10: # delete if Bayesian optimization is not finished
                    csv_file.unlink()

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {file} --format fasta --alphabet dna --max_partitions 30 --n_iter 10')
            else:
                logger.info(f'Found existing PsiPartition results for {file}')



    elif args.dataset == 'DNA-simulated':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated'
        sp_dir = Path(dna_simulated_dir) / 'PsiPartitionFast'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_simulated_dir).glob('**/*.fasta'):
            part_file = sp_dir / file.parent.name / f'{file.stem}.parts'
            part_file.parent.mkdir(exist_ok=True, parents=True)
            new_file =part_file.with_suffix(".fasta")
            os.system(f'cp {file} {new_file}')
            file = new_file

            csv_file = Path(file).with_suffix('.csv')
            if csv_file.exists():
                df = pd.read_csv(csv_file, header=None)
                if df.shape[0] != 10: # delete if Bayesian optimization is not finished
                    csv_file.unlink()

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {new_file} --format fasta --alphabet dna --max_partitions 30 --n_iter 10')
            else:
                logger.info(f'Found existing PsiPartition results for {file}')

    elif args.dataset == 'protein-empirical':
        ########################## test protein empirical ##########################
        n_iter = 20
        n_part = 20
        protein_empirical_dir = r'data/Protein/empirical'
        sp_dir = Path(protein_empirical_dir) / 'PsiPartitionFast'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(protein_empirical_dir).glob('*.fasta'):
            part_file = sp_dir / f'{file.stem}.parts'
            new_file = part_file.with_suffix(".fasta")
            # move file to PsiPartition directory
            os.system(f'cp {file} {new_file}')
            file = new_file
            
            csv_file = Path(file).with_suffix('.csv')
            if csv_file.exists():
                df = pd.read_csv(csv_file, header=None)
                if df.shape[0] != n_iter:
                    csv_file.unlink()

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {file} --format fasta --alphabet aa --max_partitions {n_part} --n_iter {n_iter}')
            else:
                logger.info(f'Found existing PsiPartition results for {file}')


    elif args.dataset == 'DNA-simulated-seqlen':
        ########################## test DNA simulated-seqlen ##########################
        dna_simulated_seqlen_dir = r'data/DNA/simulated-seqlen'
        sp_dir = Path(dna_simulated_seqlen_dir) / 'PsiPartitionFast'
        sp_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_simulated_seqlen_dir).glob('len*loci*/*.fasta'):
            part_file = sp_dir / file.parent.name / f'{file.stem}.parts'
            part_file.parent.mkdir(exist_ok=True, parents=True)
            new_file =part_file.with_suffix(".fasta")
            os.system(f'cp {file} {new_file}')
            file = new_file

            csv_file = Path(file).with_suffix('.csv')
            if csv_file.exists():
                df = pd.read_csv(csv_file, header=None)
                if df.shape[0] != 10:
                    csv_file.unlink()

            if not csv_file.exists():
                # run PsiPartition
                logger.info(f'Running PsiPartition on {file}')
                os.system(f'python {psipartition_path} --msa {new_file} --format fasta --alphabet dna --max_partitions 30 --n_iter 10')
            else:
                logger.info(f'Found existing PsiPartition results for {file}')