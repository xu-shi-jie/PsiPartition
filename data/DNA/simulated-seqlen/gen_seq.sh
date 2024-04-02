for i in `seq 1 10` ; do
    seq_len=$((1000 * i + 5000))
    n_loci=$((i+5))
    echo $seq_len
    python data/DNA/simulated-seqlen/gen_fake_dnaseq.py --n_trees 10 --n_seqs 49 --n_loci $n_loci --seq_len $seq_len --output_dir data/DNA/simulated-seqlen/
done