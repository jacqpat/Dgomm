import sys
import yaml
import random
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

################################################
################ parameters ####################
################################################

parser = op.OptionParser(conflict_handler="resolve")
parser.add_option( '-y', '--yaml', dest = 'yaml', default = '', action = 'store', help = 'Config file in yaml format needed' )
parser.add_option( '-i', '--input_name', dest = 'input_name', default='DGOMM.vcf', action = 'store', help = 'input VCF file, compressed or not.')
parser.add_option( '-o', '--output_name', dest = 'output_name', default='DGOMMd', action = 'store', help = 'Root of output files, without extension')
( options, spillover ) = parser.parse_args()
config = read_yaml(parser,options)

data_col = 9 # CHROM POS ID REF ALT QUAL INFO FORMAT
population = config['Sample_population']['pop_1']
nb_present = population[0][1]
nb_anciens = 0
for generation in population[1:]:
    nb_anciens += generation[1]
debut_adna = data_col + nb_present + 1
fin_adna = debut_adna + nb_anciens - 1

fraction_survr = config['mod_parameters']['fraction_destruction']
fraction_desam = config['mod_parameters']['fraction_deamination']
fraction_degrd = config['mod_parameters']['fraction_degradation']

vcf = options.input_name
print(vcf)
out = options.output_name + ".vcf"
print(out)

################################################
#################### main ######################
################################################

with open(vcf, "r") as fin, open(out, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        
        if random.random() < fraction_survr:
            continue
        
        columns = line.strip().split("\t")

        if columns[3] == "T":
                aa = 0
                ii = 0
                oo = 0
                for i in range(debut_adna, fin_adna):
                    if random.random() < fraction_desam:
                        columns[i] = "Y"
                for i in range(data_col+1, fin_adna):
                    if columns[i] == "Y":
                        aa += 1
                    elif columns[i] == "1":
                        ii += 1
                    else:
                        oo += 1
                if columns[4] == "C":
                    for i in range(debut_adna, fin_adna):
                        if columns[i] == "Y":
                            columns[i] = "1"
                else:
                    if aa > oo:
                        columns[4] = "C"
                        for i in range(data_col+1, fin_adna):
                            if columns[i] == "1":
                                columns[i] = "."
                            if columns[i] == "Y":
                                columns[i] = "1"
                    if aa < oo:
                        for i in range(debut_adna, fin_adna):
                            if columns[i] == "Y":
                                columns[i] = "."

        for i in range(debut_adna, fin_adna):
            if random.random() < fraction_degrd:
                columns[i] = "."

        fout.write("\t".join(columns) + "\n")
