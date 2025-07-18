import sys
import yaml
import numpy as np
import msprime
import pandas as pd
import optparse as op

################################################
################# functions ####################
################################################
def read_yaml(parser,options):
    if not options.yaml:
        parser.print_help()
        sys.exit()
    elif not options.yaml.endswith(".yaml") and not options.yaml.endswith(".yml"):
        print("File given doesn't have the YAML extension")
        sys.exit()
    yaml_file = options.yaml
    with open(yaml_file, "r") as f:
        try:
            return(yaml.safe_load(f))
        except yaml.YAMLError as e:
            print("Error during parsing :", e)

def write_custom_vcf(ts, fn):
    with open(fn, "w") as vcf:
        # Write default metadata
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        # Extract mutation sites
        for s in ts.sites():
            pos = int(s.position)
            anc = s.ancestral_state  # REF
            # Check and convert mutations mutations
            alt = set()
            for m in s.mutations:
                alt.add(m.derived_state)  # ALT
            # check there's an ALT
            if alt:
                ref = anc if anc in ["A", "T", "C", "G"] else "N"
                alt = ",".join(alt)
                # Write in VCF
                vcf.write(f"1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

################################################
################ parameters ####################
################################################

# parsing parameters (yaml and outfile)
parser = op.OptionParser(conflict_handler="resolve")
parser.add_option( '-y', '--yaml', dest = 'yaml', default = '', action = 'store', help = 'Config file in yaml format needed' )
parser.add_option( '-o', '--outfile_name', dest = 'outfile_name', default='DGOMM', action = 'store', help = 'Root of output files, without extension')
parser.add_option( '-p', '--params_file', dest = 'params_name', default='', action = 'store', help = 'TSV parameter file with three types of columns: mutation rate, rec, and population size (one for each time window)')
( options, spillover ) = parser.parse_args()
config = read_yaml(parser,options)
outfile_name = options.outfile_name
print(config)

segment_length = config['General_parameters']['L']
segment_number = config['General_parameters']['nb_seg']

# time windows
number_windows = config['Time_windows']['nb_times']
maximum_time = config['Time_windows']['t_max']
a = 0.06 # the length of time windows increases when time increases, at a speed that is determined by this coefficient

population = config['Sample_population']['pop_1']
print(population)

params=np.loadtxt(options.params_name,dtype='float',ndmin=2) # (mut,rec,pop_sizes)
print("params:",params)
nb_rep=np.shape(params)[0]
pop_init=params[0,2]
print(np.shape(params))
print("nb_rep:",nb_rep)
print("initial pop:",pop_init)

# computation of time windows based on the above parameters
times = -np.ones(shape = number_windows, dtype='float')
for i in range(number_windows):
     times[i]=(np.exp(np.log(1+a*maximum_time)*i/(number_windows-1))-1)/a
print ("Population size changes at the following times (in generations):")
print (times)
print ("")


# computation of samples list
pop = []
for i in range(len(population)):
    gen = population[i][0]
    nbr = population[i][1]
    pop.append(msprime.SampleSet(nbr, time = gen, ploidy = 1))
print(pop)

################################################
#################### main ######################
################################################

# generate samples
for rep in range(nb_rep):
    rep_params=params[rep,:]
    # creates demographic events
    mymut=rep_params[0]
    myrec=rep_params[1]
    N=rep_params[2]
    mydemo=[]
    demography = msprime.Demography()
    demography.add_population(initial_size=pop_init)
    for i in range(1,len(times)):
        demography.add_population_parameters_change(time=times[i],initial_size= rep_params[2+i])
    # simulate each segment   
    for i in range(segment_number):
        #tree_sequence = msprime.simulate(sample_size=n, Ne=N, length=L, recombination_rate=myrec, mutation_rate=mymut, demographic_events=mydemo)
        tree_sequence = msprime.sim_ancestry(samples=pop,
                            recombination_rate=myrec,
                            sequence_length=segment_length,
                            demography=demography
                        )
    # add mutations
    mutated_tree_sequence = msprime.sim_mutations(tree_sequence, rate=mymut)
	# export (to be adapted to get the format required by SMC, specific functions probably exist)
    with open(outfile_name+".vcf", "w") as vcf_file:
        mutated_tree_sequence.write_vcf(vcf_file)
