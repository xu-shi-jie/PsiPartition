import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

aa_str = 'ARNDCQEGHILKMFPSTWYV-'
dna_str = 'ACGT-'


def extract_stat(p):
    """ Extract statistics from IQ-TREE iqtree file.
    Args:
        p (str): Path to the iqtree file.
    Returns:
        stat: Statistics inlcuding log-likelihood, AIC, AICc, BIC, and number of parameters.
    """
    content = open(p).read()
    ll = re.search(r'Log-likelihood of the tree:\s+(-?\d+\.\d+)', content)
    ll = float(ll.group(1))
    aic = re.search(
        r'Akaike information criterion \(AIC\) score:\s+(-?\d+\.\d+)', content)
    aic = float(aic.group(1))
    aicc = re.search(
        r'Corrected Akaike information criterion \(AICc\) score:\s+(-?\d+\.\d+)', content)
    aicc = float(aicc.group(1))
    bic = re.search(
        r'Bayesian information criterion \(BIC\) score:\s+(-?\d+\.\d+)', content)
    bic = float(bic.group(1))
    num_params = re.search(
        r'Number of free parameters \(#branches \+ #model parameters\):\s+(\d+)', content)
    num_params = int(num_params.group(1))
    return {
        'log-likelihood': ll,
        'AIC': aic,
        'AICc': aicc,
        'BIC': bic,
        'num_params': num_params,
    }


def default_iqtree(options, log_file):
    """ Run IQ-TREE with options, and write log to log_file.
        Only 'iqtree', '-redo', '-safe' are included by default.

    Args:
        options (dict): Dictionary of options.
        log_file (file handle): file handle to write log to.
    """
    # use -te *.treefile to specify a starting tree and do NOT optimize?
    # use -t *.treefile to specify a starting tree and optimize
    # use -wsl to compute the site-wise log-likelihood
    # use -m to specify a model
    subprocess.call(['./ThirdParty/iqtree-1.6.12-Linux/bin/iqtree', '-redo', '-safe',
                    *options], stderr=log_file, stdout=log_file,)


def remove_caches(file_path):
    """ Remove IQ-TREE caches.

    Args:
        file_path (str): Path to the iqtree file.
    """
    for file in Path(file_path).parent.glob(file_path.stem + '.*'):
        if file.suffix in [
            '.gz', '.log', '.mldist', '.bionj', '.best_scheme',
                '.nex', '.treefile', '.sitelh', '.rates']:
            file.unlink()

def write_part_file(part_file, indices):
    """ Write partition file.

    Args:
        part_file (str): Path to the partition file.
        indices (list): List of indices.
    """
    with open(part_file, 'w') as f:
        f.write('#nexus\nbegin sets;\n')
        parts = np.unique(indices)
        for p in parts:
            f.write(
                f'charset part{p} = {" ".join([str(i+1) for i in np.where(indices == p)[0]])};\n')
        f.write('end;')


def convert_RP_to_nexus(infile, outfile):
    content = open(infile).read()
    phylip_rp = re.search(
        r'PHYLIP  style(.+)', content, re.DOTALL | re.MULTILINE).group(1).strip()
    with open(outfile, 'w') as f:
        f.write('#nexus\nbegin sets;\n')
        for line in phylip_rp.split('\n'):
            line = line.strip()
            if line == '':
                continue
            parts = line.split()
            f.write('charset ' + parts[1] + ''.join(parts[2:]).replace(',', ' ') + ';\n')
        f.write('end;')


def is_one_state_error(infile) -> bool:
    content = open(infile).read()
    if 'ERROR: Only one state is observed' in content:
        return True
    return False

def transfer_distance(s1, s2):
    """ Transfer distance of two list s1 and s2."""
    # we consider s1 and s2 are the partitions of the same set.
    # cost matrix for the Hungarian algorithm
    cost_matrix = np.zeros((len(s1), len(s2)))
    for i, x in enumerate(s1):
        for j, y in enumerate(s2):
            cost_matrix[i, j] = len(x.intersection(y))
    row_ind, col_ind = linear_sum_assignment(-cost_matrix)
    return 1 - cost_matrix[row_ind, col_ind].sum() / sum(map(len, s1))

if __name__ == '__main__':
    s1 = [{1,2,3}, {4,5,6}]
    s2 = [{1,6}, {2,5,}, {3,4}]
    print(transfer_distance(s1, s2))