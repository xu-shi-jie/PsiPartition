import argparse
import os
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from binarytree import Node
from Bio import AlignIO
from tqdm import tqdm


def gen_bin_tree(n):
    nodes = [Node(i) for i in range(n)]
    i = n
    while len(nodes) > 1:
        random.shuffle(nodes)
        # merge the last two nodes recursively
        new_node = Node(i, nodes.pop(), nodes.pop())
        nodes.append(new_node)
        i += 1
    return new_node


def gen_newick(node):
    def f():
        return np.random.default_rng().exponential(scale=1)
    if node is None:
        return ''
    if node.left is None and node.right is None:
        return 'node' + str(node.value)
    return f'({gen_newick(node.left)}:{f()},{gen_newick(node.right)}:{f()})'


def gen_seq(args, tree_path, output_path):
    freqs = np.random.random(4)
    freqs = freqs / np.sum(freqs)

    os.system(
        './' + args.seqgen_path + \
        ' -r'+' '.join([f'{x:.3f}' for x in np.random.random(6)]) + \
        ' -f'+' '.join([f'{x:.3f}' for x in freqs]) + \
        ' -m GTR -g 4 -l ' + str(args.seq_len // args.n_loci) + \
        ' -n 1 < ' + tree_path + ' > ' + output_path

    )

def concat_fasta(fasta_paths, output_path):
    seqs = defaultdict(str)
    for f in fasta_paths:
        for r in AlignIO.read(f, 'phylip'):
            seqs[r.id] += str(r.seq)
        
    with open(output_path, 'w') as f1:
        # write to fasta
        for k, v in seqs.items():
            f1.write(f'>{k}\n{v}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_trees', type=int, default=10)
    parser.add_argument('--n_seqs', type=int, default=49)
    parser.add_argument('--n_loci', type=int, default=10)
    parser.add_argument('--seq_len', type=int, default=5000)
    parser.add_argument('--output_dir', type=str, default='data/DNA/simulated')
    parser.add_argument('--seqgen_path', type=str,
                        default='ThirdParty/Seq-Gen/source/seq-gen')
    args = parser.parse_args()

    # shutil.rmtree(args.output_dir, ignore_errors=True)
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    args.output_dir = args.output_dir + '/' + f'len_{args.seq_len}_loci_{args.n_loci}'
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    for i in tqdm(range(args.n_trees)):
        random_id = str(i)
        tree = gen_bin_tree(args.n_seqs)
        tree_path = f'{args.output_dir}/{random_id}.tree'
        with open(tree_path, 'w') as f1:
            newick = gen_newick(tree) + ';'
            f1.write(newick)
        outputs = []
        for j in range(args.n_loci):
            output_path = f'{args.output_dir}/{random_id}_{j}.phylip'
            outputs.append(output_path)
            gen_seq(args, tree_path, output_path)
        concat_fasta(outputs, f'{args.output_dir}/{random_id}.fasta')
        for f in outputs:
            os.remove(f)