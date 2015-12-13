# About VBBMM #

VBBMM performs variational bayes binomial mixture model clustering on genome sequencing allelic count data in n-dimensions.

# How to Install #

VBBMM is a python script and thus requires no installation. However, it does depend on [PyVMP](https://bitbucket.org/aroth85/pyvmp/wiki/Home). Your python will need to have the PyVMP library installed in order to run VBBMM.

# Input File

VBBMM takes a tsv input file containing the following required columns:

* chr: Chromosome. This is not used in the cluster, but is used as a identifier.
* pos: Position. This is not used in the cluster, but is used as a identifier.
* T1_TR: Number of reads supporting the tumor reference allele in the first dimension
* T1_TA: Number of reads supporting the tumor variant allele in the first dimension

To add additional dimensions to cluster on, you simply have to add at Ti_TR and Ti_TA for each new ith dimension (e.g. T2_TR, T2_TA)

# How to Run

To run vbbmm, for instance on 2 dimensions, you can use the following command:

```
#!bash

python vbbmm.py --in_file in.tsv --out_file out.tsv --num_samples 2
```


