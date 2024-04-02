for file in $(ls data/BayesOpt/*.fasta); do
    echo $file
    python comparison/test_bayes_opt.py --msa $file --n_iter 60 --opt bayes
    python comparison/test_bayes_opt.py --msa $file --n_iter 1000 --opt random 
done