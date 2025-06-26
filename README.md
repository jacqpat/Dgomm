# Dgomm
**DGOMM** (**D**egradation **G**enerator **O**f **M**utated **M**Sprime) is a simple Python pipeline meant to simulate post-mortem degradation of DNA, such as fragmentation and deamination (_T -> C_), in a basic and easy to use manner. 

DGOMM' use as input a **VCF** (_Variant Calling Format_), generated from real or simulated data. The pipeline comes with an optional script to simulate such a file using a pre-set demographic history and the population genetics simulator **MSprime** (https://tskit.dev/msprime/docs/stable/intro.html).

Once you have a VCF, we can set the degradation and filtration parameters. "Degradation parameters" here refer to the necessary variables for the simulation of post-mortem DNA mutation and degradation's consequences. In this first version of the pipeline, those tools are very basic: for each of the VCF's entries, an ancient DNA sample has an odd of losing the information at this position. If the entry's reference allele is T, then Deamination is applied to randomly transform some of the ancient entries into C.

## The Pipeline
### (Optional 1) Set base demographic parameters
### (Optional 2) Generate VCF with MSprime
### 1) Set degradation and filtration parameters
### 2) Degrade VCF
### 3) Filter degraded VCF
