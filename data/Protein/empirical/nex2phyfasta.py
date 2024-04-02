from pathlib import Path

from Bio import AlignIO

for file in Path('data/Protein/empirical').glob('*.nex'):
    aln = AlignIO.read(file, 'nexus')
    print(f'Found {aln.get_alignment_length()} sites and {len(aln)} sequences in {file}.')
    # replace spaces in sequence names with underscores
    for seq in aln:
        seq.id = seq.id.replace(' ', '_')
    # convert to fasta and save
    AlignIO.write(aln, file.with_suffix('.fasta'), 'fasta')
    with open(file.with_suffix('.phy'), 'w') as f:
        f.write(f'{len(aln)}\t{aln.get_alignment_length()}\n')
        for r in aln:
            f.write(f'{r.id}\t{str(r.seq)}\n')