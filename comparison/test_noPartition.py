import os
import sys
from pathlib import Path

sys.path.append('.')
import argparse

from loguru import logger

from utils import default_iqtree, remove_caches

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=str, default='protein-simulated')
    args = parser.parse_args()
    
    if args.dataset == 'DNA-empirical':
        ########################## test DNA empirical ##########################
        dna_empirical_dir = r'data/DNA/empirical'
        np_dir = Path(dna_empirical_dir) / 'noPartition'
        np_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_empirical_dir).glob('*.fasta'):
            logger.info(f'Running IQ-TREE on {file}')
            remove_caches(np_dir / f'{file.stem}.fasta')
            default_iqtree(['-s', file, '-nt', '4', '-pre', np_dir / file.stem], open(np_dir / f'{file.stem}.log', 'w'))
    elif args.dataset == 'DNA-simulated':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated'
        np_dir = Path(dna_simulated_dir) / 'noPartition'
        np_dir.mkdir(exist_ok=True, parents=True)
        for file in Path(dna_simulated_dir).glob('loci*/*.fasta'):
            (np_dir / file.parent.name).mkdir(exist_ok=True, parents=True)
            logger.info(f'Running IQ-TREE on {file}.')
            remove_caches(np_dir / f'{file.stem}.fasta')
            default_iqtree(['-s', file, '-nt', '4', '-pre', np_dir / file.parent.name / file.stem, '-bb', '1000'], open(np_dir / f'{file.stem}.log', 'w'))
    elif args.dataset == 'protein-empirical':
        ########################## test protein empirical ##########################
        protein_empirical_dir = r'data/Protein/empirical'
        np_dir = Path(protein_empirical_dir) / 'noPartition'
        np_dir.mkdir(exist_ok=True, parents=True)
        for file in Path(protein_empirical_dir).glob('*.fasta'):
            logger.info(f'Running IQ-TREE on {file}.')
            remove_caches(np_dir / f'{file.stem}.fasta')
            default_iqtree(['-s', file, '-nt', '4', '-pre', np_dir / file.stem], open(np_dir / f'{file.stem}.log', 'w'))
    elif args.dataset == 'protein-simulated':
        ########################## test protein simulated ##########################
        protein_simulated_dir = r'data/Protein/simulated'
        np_dir = Path(protein_simulated_dir) / 'noPartition'
        np_dir.mkdir(exist_ok=True, parents=True)
        for file in Path(protein_simulated_dir).glob('loci*/*.fasta'):
            (np_dir / file.parent.name).mkdir(exist_ok=True, parents=True)
            logger.info(f'Running IQ-TREE on {file}.')
            remove_caches(np_dir / f'{file.stem}.fasta')
            default_iqtree(['-s', file, '-nt', '4', '-pre', np_dir / file.parent.name / file.stem, '-bb', '1000'], open(np_dir / f'{file.stem}.log', 'w'))
    elif args.dataset == 'DNA-simulated-seqlen':
        ########################## test DNA-seqlen simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated-seqlen'
        np_dir = Path(dna_simulated_dir) / 'noPartition'
        np_dir.mkdir(exist_ok=True, parents=True)
        for file in Path(dna_simulated_dir).glob('len*loci*/*.fasta'):
            (np_dir / file.parent.name).mkdir(exist_ok=True, parents=True)
            logger.info(f'Running IQ-TREE on {file}.')
            remove_caches(np_dir / f'{file.stem}.fasta')
            default_iqtree(['-s', file, '-nt', '8', '-pre', np_dir / file.parent.name / file.stem, ], open(np_dir / f'{file.stem}.log', 'w'))