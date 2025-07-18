# Dgomm
**DGOMM** (**D**egradation **G**enerator **O**f **M**utated **M**Sprime) is a simple Python pipeline meant to simulate post-mortem degradation of DNA, such as fragmentation and deamination (_T -> C_), in a basic and easy to use manner. 

DGOMM' use as input a **VCF** (_Variant Calling Format_), generated from real or simulated data. The pipeline comes with an optional script to simulate such a file using a pre-set demographic history and the population genetics simulator **MSprime** (https://tskit.dev/msprime/docs/stable/intro.html).

Once you have a VCF, we can set the degradation and filtration parameters. "Degradation parameters" here refer to the necessary variables for the simulation of post-mortem DNA mutation and degradation's consequences. In this first version of the pipeline, those tools are very basic: for each of the VCF's entries, an ancient DNA sample has an odd of losing the information at this position. If the entry's reference allele is T, then Deamination is applied to randomly transform some of the ancient entries into C.
## Installation
DGOMM is a very simple collection of python scripts. Simply download the zipped files, decompress, and use command lines from a shell to get the degraded VCF you want.
## Dependencies
Python 3 pipeline
Required modules:
- random (standard library)
- optparse (standard library)
## The Pipeline
### (Optional 1) Set base demographic parameters
If you do not have a VCF file yet and wish to get one that match a specific population scenario for a population, you must first create a TSV and YAML files.
The TSV file (preferably with the extension ".params") detail the size of your theoretical population and its growth (or contraction) throughout the generations. It is very simple, as you can see in the given example const_small.params:

| 1e-8 | 5e-9 | 1000 | 1000 | 1000 | [...] | 1000 |

With first the **Mutation Rate**, then the **Recombination Rate**, and finally as many columns as there will be **Time Windows** in your simulation. For each time window, the value given is the **Ne**, the Effective Population, for that era. In the example .params given, the simulated population was always constant with 1000 effective reproducers no matter the era. But you are free to make your own .params file with your own scenarios (constant growth, reduction, constant population then bottleneck...) The third column is the first, and therefore closest to modern, time window.
The YAML file detail the various parameters needed for MSprime to simulate your population. For now, the paramaters included are:
In General:
- **L**: the Length (in pb) of the sequences to simulate.
- **nb_seg**: the number of sequences (segments) to simulate.
For the time windows:
- **nb_times**: number of time windows. Must match the number of column (minus two) in the params file.
- **t_max**: the age of the oldest time window.
- **a**: the length of time windows increases when time increases, and this coefficient determine the speed of the process.
And finally for the sample population:
- **pop_1**: the number of individuals "sampled" for each generation. The format is such: [0,15] with the first number being the generation (here "0" ie the present) and the second the number of samples. Normally, it should work no matter the order of the generations within the list, but it has not been tested so to be safe keep it so there's the most recent generation in first and the oldest in last.

### (Optional 2) Generate VCF with MSprime

Once you've create both the .params and .yml files, you can pass them as parameters to **simul_genome.py** like so:

python simul_genome.py -y example.yml -p example.params -o examplar

output parameter ("-o") is for naming the VCF output and is optional but recommended.

### 1) Set degradation and filtration parameters
### 2) Degrade VCF
### 3) Filter degraded VCF
