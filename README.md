# PsiPartition
![](logo.png)

This is the official implementation of the paper "PsiPartition: Improved Site Partitioning for Genomic Data by Parameterized Sorting Indices and Bayesian Optimization" [Biorxiv](https://www.biorxiv.org/cgi/content/short/2024.04.03.588030v1).

## Installation

```bash
conda create -n PsiPartition
conda activate PsiPartition
pip install -r requirements.txt
```

To enable the Bayesian optimization, you need to register an account on [Weights & Biases](https://wandb.ai/).
## Usage

```bash
python PsiPartition_wandb.py --msa <MSA File> --format <fasta or phylip> --alphabet <dna or aa> --max_partitions <max_partitions> --n_iter <number of iterations>
```

## Help
If you have any questions, please contact me at `shijie.xu@ees.hokudai.ac.jp`.
