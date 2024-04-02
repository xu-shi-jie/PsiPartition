import os
import sys
from pathlib import Path

sys.path.append('.')
import argparse

from loguru import logger

from utils import (convert_RP_to_nexus, default_iqtree, is_one_state_error,
                   remove_caches)

if __name__ == '__main__':
    tiger_path = r'./ThirdParty/mPartition/tiger_original/tiger'


    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=str, default='DNA-simulated')
    parser.add_argument('--div_factor', type=int, default=4)
    args = parser.parse_args()
    
    div_factor = args.div_factor
    assert isinstance(div_factor, int)
    suffix = f'_{div_factor}.0.txt'

    if args.dataset == 'DNA-empirical':
        ########################## test DNA empirical ##########################
        dna_empirical_dir = r'data/DNA/empirical'
        rate_dir = Path(dna_empirical_dir) / f'RatePartition-{div_factor}'
        rate_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_empirical_dir).glob('*.fasta'):
            rate_file = rate_dir / f'{file.stem}.rates'
            if rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Found existing IQ-TREE results for {rate_file}')
                continue

            if not rate_file.exists():
                logger.info(f'Running TIGER on {file}')
                with open('tmp_RatePartition_dna_empirical.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'{tiger_path} -in {file} -rl {rate_file}\n')
                os.system('bash -i ./tmp_RatePartition_dna_empirical.sh')
                os.system('rm ./tmp_RatePartition_dna_empirical.sh')
            
            partition_file = Path(str(rate_file) + suffix)
            if not partition_file.exists(): 
                with open('tmp_RatePartition_dna_empirical.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'python ThirdParty/RatePartitions/rate_partitions.py {rate_file} {div_factor}\n')
                os.system('bash -i ./tmp_RatePartition_dna_empirical.sh')
                os.system('rm ./tmp_RatePartition_dna_empirical.sh')

            if not rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Running IQ-TREE on {file}')
                remove_caches(rate_file)
                convert_RP_to_nexus(partition_file, partition_file.with_suffix('.nex'))
                default_iqtree(['-s', file, '-spp', partition_file.with_suffix('.nex'), '-nt', '4', '-pre', rate_file.with_suffix('')], open(rate_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {rate_file}')

    elif args.dataset == 'DNA-simulated':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated'
        rate_dir = Path(dna_simulated_dir) / f'RatePartition-{div_factor}'
        rate_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_simulated_dir).glob('loci*/*.fasta'):
            rate_file = rate_dir / file.parent.name / f'{file.stem}.rates'
            (rate_file.parent).mkdir(exist_ok=True, parents=True)

            if rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Found existing IQ-TREE results for {rate_file}')
                continue

            if not rate_file.exists():
                logger.info(f'Running TIGER on {file}')
                with open(f'tmp_RatePartition_dna_simulated-{div_factor}.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'{tiger_path} -in {file} -rl {rate_file}\n')
                os.system(f'bash -i ./tmp_RatePartition_dna_simulated-{div_factor}.sh')
                os.system(f'rm ./tmp_RatePartition_dna_simulated-{div_factor}.sh')

            partition_file = Path(str(rate_file) + suffix)
            if not partition_file.exists(): 
                with open(f'tmp_RatePartition_dna_simulated-{div_factor}.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'python ThirdParty/RatePartitions/rate_partitions.py {rate_file} {div_factor}\n')
                os.system(f'bash -i ./tmp_RatePartition_dna_simulated-{div_factor}.sh')
                os.system(f'rm ./tmp_RatePartition_dna_simulated-{div_factor}.sh')

            if not rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Running IQ-TREE on {file}')
                remove_caches(rate_file)
                convert_RP_to_nexus(partition_file, partition_file.with_suffix('.nex'))
                default_iqtree(['-s', file, '-spp', partition_file.with_suffix('.nex'), '-nt', '4', '-pre', rate_file.with_suffix('')], open(rate_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {rate_file}')

    elif args.dataset == 'protein-empirical':
        ########################## test protein empirical ##########################
        protein_empirical_dir = r'data/Protein/empirical'
        rate_dir = Path(protein_empirical_dir) / f'RatePartition-{div_factor}'
        rate_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(protein_empirical_dir).glob('*.fasta'):
            rate_file = rate_dir / f'{file.stem}.rates'
            if rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Found existing IQ-TREE results for {rate_file}')
                continue

            if not rate_file.exists():
                logger.info(f'Running TIGER on {file}')
                with open('tmp_RatePartition_protein_empirical.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'{tiger_path} -in {file} -rl {rate_file}\n')
                os.system('bash -i ./tmp_RatePartition_protein_empirical.sh')
                os.system('rm ./tmp_RatePartition_protein_empirical.sh')
            
            partition_file = Path(str(rate_file) + suffix)
            if not partition_file.exists(): 
                with open('tmp_RatePartition_protein_empirical.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'python ThirdParty/RatePartitions/rate_partitions.py {rate_file} {div_factor}\n')
                os.system('bash -i ./tmp_RatePartition_protein_empirical.sh')
                os.system('rm ./tmp_RatePartition_protein_empirical.sh')

            if not rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Running IQ-TREE on {file}')
                remove_caches(rate_file)
                convert_RP_to_nexus(partition_file, partition_file.with_suffix('.nex'))
                default_iqtree(['-s', file, '-spp', partition_file.with_suffix('.nex'), '-nt', '4', '-pre', rate_file.with_suffix('')], open(rate_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {rate_file}')

    elif args.dataset == 'DNA-simulated-seqlen':
        ########################## test DNA-seqlen simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated-seqlen'
        rate_dir = Path(dna_simulated_dir) / f'RatePartition-{div_factor}'
        rate_dir.mkdir(exist_ok=True, parents=True)

        for file in Path(dna_simulated_dir).glob('len*loci*/*.fasta'):
            rate_file = rate_dir / file.parent.name / f'{file.stem}.rates'
            (rate_file.parent).mkdir(exist_ok=True, parents=True)

            if rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Found existing IQ-TREE results for {rate_file}')
                continue

            if not rate_file.exists():
                logger.info(f'Running TIGER on {file}')
                with open(f'tmp_RatePartition_dna_simulated-seqlen-{div_factor}.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'{tiger_path} -in {file} -rl {rate_file}\n')
                os.system(f'bash -i ./tmp_RatePartition_dna_simulated-seqlen-{div_factor}.sh')
                os.system(f'rm ./tmp_RatePartition_dna_simulated-seqlen-{div_factor}.sh')

            partition_file = Path(str(rate_file) + suffix)
            if not partition_file.exists(): 
                with open(f'tmp_RatePartition_dna_simulated-seqlen-{div_factor}.sh', 'w') as f:
                    f.write('conda activate py2.7\n')
                    f.write(f'python ThirdParty/RatePartitions/rate_partitions.py {rate_file} {div_factor}\n')
                os.system(f'bash -i ./tmp_RatePartition_dna_simulated-seqlen-{div_factor}.sh')
                os.system(f'rm ./tmp_RatePartition_dna_simulated-seqlen-{div_factor}.sh')

            if not rate_file.with_suffix('.iqtree').exists():
                logger.info(f'Running IQ-TREE on {file}')
                remove_caches(rate_file)
                convert_RP_to_nexus(partition_file, partition_file.with_suffix('.nex'))
                default_iqtree(['-s', file, '-spp', partition_file.with_suffix('.nex'), '-nt', '8', '-pre', rate_file.with_suffix('')], open(rate_file.with_suffix('.log'), 'w'))