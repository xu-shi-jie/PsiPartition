from pathlib import Path

from Bio import AlignIO

# convert fasta to phylip
for file in Path('data/DNA/simulated-seqlen').glob('len*loci*/*.fasta'):
    rec = AlignIO.read(file, 'fasta')
    with open(file.with_suffix('.phy'), 'w') as f:
        n_seqs, n_sites = len(rec), rec.get_alignment_length()
        f.write(f'{n_seqs}\t{n_sites}\n')
        for seq in rec:
            f.write(f'{seq.id}\t{seq.seq}\n')