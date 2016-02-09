# About VBBMM #

VBBMM performs variational bayes binomial mixture model clustering on genome sequencing allelic count data in n-dimensions.

# How to Install #

VBBMM is a python script and thus requires no installation. However, it does depend on [PyVMP](https://bitbucket.org/aroth85/pyvmp/wiki/Home). Your python will need to have the PyVMP library installed in order to run VBBMM.

# Input File

VBBMM takes a tsv input file containing the following required columns:

1. chr: Chromosome. This is not used in the cluster, but is used as a identifier.
1. pos: Position. This is not used in the cluster, but is used as a identifier.
1. T1_TR: Number of reads supporting the reference allele in the first dimension
1. T1_TA: Number of reads supporting the variant allele in the first dimension

To add additional dimensions to cluster on, you simply have to add a Ti_TR and a Ti_TA column for each new ith dimension (e.g. T2_TR, T2_TA)

# How to Run

To run VBBMM, for instance on 2 dimensions, you can use the following command:

```bash
python vbbmm.py --in_file in.tsv --out_file out.tsv --num_samples 2
```

Other optional parameters to consider are:

* --num_components: Number of components (clusters) to consider
* --mix_weight_priors: Dirichlet prior on component mix-weights. Default is 1e-6.
* --min_depth: Minimum depth in any sample required to include position for analysis.
* --max_iters: Maximum number of VB iterations to do. Default is 100.

# Output

The output of VBBMM will be a tsv file with:

1. chr: Identifier from input tsv file. 
1. pos: Identifier from input tsv file.
1. cluster: Cluster identifier which this data point belongs to.
