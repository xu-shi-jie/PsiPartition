# mPartition requires python 2.7
# do NOT run this script unless you know what you are doing

import os
import sys

sys.path.append('.')
import argparse
from pathlib import Path

from loguru import logger

from utils import default_iqtree, remove_caches

if __name__ == '__main__':
    mp_path = r'ThirdParty/mPartition/mPartition.py'

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=str, default='protein-empirical')
    args = parser.parse_args()
    
    if args.dataset == 'DNA-empirical':
        ########################## test DNA empirical ##########################
        dna_empirical_dir = r'data/DNA/empirical'
        (Path(dna_empirical_dir) / 'mPartition').mkdir(exist_ok=True, parents=True)
        for file in Path(dna_empirical_dir).glob('*.phy'):
            par_file = Path('ThirdParty/mPartition/Results') / f'par.{file.stem}.phy'
            new_par_file = Path(dna_empirical_dir) / 'mPartition' / par_file.name
            if par_file.exists():
                logger.info(f'Found existing mPartition results for {par_file}')
                par_file.rename(new_par_file)
            elif new_par_file.exists():
                logger.info(f'Found existing mPartition results for {new_par_file}')
            else:
                logger.info(f'Testing mPartition on {file}')
                with open('tmp_mPartition_dna_empirical.sh', 'w') as f:
                    f.write('cd ThirdParty/mPartition\n')
                    f.write('conda activate py2.7\n')
                    f.write(f'python mPartition.py -f ../../{file}\n')
                    f.write(f'mv Results/par.{file.stem}.phy ../../{dna_empirical_dir}/mPartition\n')
                os.system('bash -i ./tmp_mPartition_dna_empirical.sh')
                os.system('rm ./tmp_mPartition_dna_empirical.sh')
            
            par_file = new_par_file
            iqtree_file = par_file.with_suffix('.iqtree')
            if not iqtree_file.exists():
                logger.info(f'Testing IQ-TREE on {par_file}')
                remove_caches(par_file)
                default_iqtree(['-s', file, '-spp', par_file, '-nt', '4', '-pre', par_file.with_suffix('')], open(par_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {par_file}')

    elif args.dataset == 'DNA-simulated':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated'
        (Path(dna_simulated_dir) / 'mPartition').mkdir(exist_ok=True, parents=True)
        for file in Path(dna_simulated_dir).glob('loci*/*.phy'):
            (Path(dna_simulated_dir) / 'mPartition' / file.parent.name).mkdir(exist_ok=True, parents=True)
            par_file = Path('ThirdParty/mPartition/Results') / f'par.{file.stem}.phy'
            new_par_file = Path(dna_simulated_dir) / 'mPartition' / file.parent.name / par_file.name
            if par_file.exists():
                logger.info(f'Found existing mPartition results for {par_file}')
                par_file.rename(new_par_file)
            elif new_par_file.exists():
                logger.info(f'Found existing mPartition results for {new_par_file}')
            else:
                logger.info(f'Testing mPartition on {file}')
                with open('tmp_mPartition_dna_simulated.sh', 'w') as f:
                    f.write('cd ThirdParty/mPartition\n')
                    f.write('conda activate py2.7\n')
                    f.write(f'python mPartition.py -f ../../{file}\n')
                    f.write(f'mv Results/par.{file.stem}.phy ../../{dna_simulated_dir}/mPartition/{file.parent.name}\n')
                os.system('bash -i ./tmp_mPartition_dna_simulated.sh')
                os.system('rm ./tmp_mPartition_dna_simulated.sh')

            par_file = new_par_file
            iqtree_file = par_file.with_suffix('.iqtree')
            if not iqtree_file.exists():
                logger.info(f'Testing IQ-TREE on {par_file}')
                remove_caches(par_file)
                default_iqtree(['-s', file, '-spp', par_file, '-nt', '4', '-pre', par_file.with_suffix('')], open(par_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {par_file}')


    elif args.dataset == 'protein-empirical':
        ########################## test protein empirical ##########################
        protein_empirical_dir = r'data/Protein/empirical'
        (Path(protein_empirical_dir) / 'mPartition').mkdir(exist_ok=True, parents=True)
        for file in Path(protein_empirical_dir).glob('*.phy'):
            par_file = Path('ThirdParty/mPartition/Results') / f'par.{file.stem}.phy'
            new_par_file = Path(protein_empirical_dir) / 'mPartition' / par_file.name
            if par_file.exists():
                logger.info(f'Found existing mPartition results for {par_file}')
                par_file.rename(new_par_file)
            elif new_par_file.exists():
                logger.info(f'Found existing mPartition results for {new_par_file}')
            else:
                logger.info(f'Testing mPartition on {file}')
                with open('tmp_mPartition_protein_empirical.sh', 'w') as f:
                    f.write('cd ThirdParty/mPartition\n')
                    f.write('conda activate py2.7\n')
                    f.write(f'python mPartition.py -f ../../{file}\n')
                    f.write(f'mv Results/par.{file.stem}.phy ../../{protein_empirical_dir}/mPartition\n')
                os.system('bash -i ./tmp_mPartition_protein_empirical.sh')
                os.system('rm ./tmp_mPartition_protein_empirical.sh')
            
            par_file = new_par_file
            iqtree_file = par_file.with_suffix('.iqtree')
            if not iqtree_file.exists():
                logger.info(f'Testing IQ-TREE on {par_file}')
                remove_caches(par_file)
                default_iqtree(['-s', file, '-spp', par_file, '-nt', '4', '-pre', par_file.with_suffix('')], open(par_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {par_file}')

    elif args.dataset == 'DNA-simulated-seqlen':
        ########################## test DNA simulated ##########################
        dna_simulated_dir = r'data/DNA/simulated-seqlen'
        (Path(dna_simulated_dir) / 'mPartition').mkdir(exist_ok=True, parents=True)
        for file in Path(dna_simulated_dir).glob('len*loci*/*.phy'):
            (Path(dna_simulated_dir) / 'mPartition' / file.parent.name).mkdir(exist_ok=True, parents=True)
            par_file = Path('ThirdParty/mPartition/Results') / f'par.{file.stem}.phy'
            new_par_file = Path(dna_simulated_dir) / 'mPartition' / file.parent.name / par_file.name
            if par_file.exists():
                logger.info(f'Found existing mPartition results for {par_file}')
                par_file.rename(new_par_file)
            elif new_par_file.exists():
                logger.info(f'Found existing mPartition results for {new_par_file}')
            else:
                logger.info(f'Testing mPartition on {file}')
                with open('tmp_mPartition_dna_simulated.sh', 'w') as f:
                    f.write('cd ThirdParty/mPartition\n')
                    f.write('conda activate py2.7\n')
                    f.write(f'python mPartition.py -f ../../{file}\n')
                    f.write(f'mv Results/par.{file.stem}.phy ../../{dna_simulated_dir}/mPartition/{file.parent.name}\n')
                os.system('bash -i ./tmp_mPartition_dna_simulated.sh')
                os.system('rm ./tmp_mPartition_dna_simulated.sh')

            par_file = new_par_file
            iqtree_file = par_file.with_suffix('.iqtree')
            if not iqtree_file.exists():
                logger.info(f'Testing IQ-TREE on {par_file}')
                remove_caches(par_file)
                default_iqtree(['-s', file, '-spp', par_file, '-nt', '4', '-pre', par_file.with_suffix('')], open(par_file.with_suffix('.log'), 'w'))
            else:
                logger.info(f'Found existing IQ-TREE results for {par_file}')