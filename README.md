# About VBBMM #

VBBMM is a python script that performs variational bayes binomial mixture model clustering on allele count data in n-dimensions.

# Dependencies #

VBBMM depends on [PyVMP](https://bitbucket.org/aroth85/pyvmp/wiki/Home).

# How to Use

VBBMM takes as tsv input file containing the following required columns:

* chr: Chromosome
* pos: Position
* Ti_TR: Number of reads supporting the tumor reference allele in the ith dimension
* Ti_TA: Number of reads supporting the tumor variant allele in the ith dimension

