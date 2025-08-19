import os, sys
import subprocess
from multiprocessing import Pool
from dotenv import dotenv_values
import pandas as pd
import argparse

from determine_ibd import run_transform as ibd

config = dotenv_values('../.env') # This will change when the repos are merged 

##################################################

def execute_command(command):
    with open(os.devnull, 'w') as devnull:
        subprocess.run(command, shell=True, stdout=devnull, stderr=devnull)

def execute_commands(commands, n_processes):
    with Pool(n_processes) as pool:
        pool.map(execute_command, commands)

##################################################

# all stuff from the .env file

data_dir = config['DATA_DIR']
gmap_dir = config['GMAP_DIR'] # genetic map files 
shapeit5 = config['SHAPEIT5_EXE'] # SHAPEIT5 executable
germline2 = config['GERMLINE2_EXE'] # GERMLINE2 executable
plink2 = config['PLINK2_EXE'] # PLINK2 executable
python3 = config['PYTHON3_EXE'] # Python3 executable 

################################################## 


'''

EXAMPLE: 

python3 simulate_shared_segments.py --vcf uniform3_size20_sim999_all_chr_qced.vcf.gz

'''

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', type=str, required=True, help='Input VCF file')
parser.add_argument('-s', '--size', type=int, default=20, help='Total number of simulated individuals')
parser.add_argument('-p', '--population', type=str, default='EUR', help='1KG superpopulation')
parser.add_argument('-m', '--missing', type=str, help='Comma-separated missing IDs to exclude from simulation')
args = parser.parse_args()


#vcf_file = sys.argv[1] 
vcf_file = args.vcf
#print (vcf_file)
working_dir_bits = vcf_file.split('/')
filehandle = working_dir_bits[-1].split('.')[0]
filehandlefull = working_dir_bits[-1]
working_dir_bits = working_dir_bits[:-1]
working_dir = '/'.join(working_dir_bits)

### NEW STEP: MAKE AN INSIDE DIR for ss stuff and attach to working_dir 
try:
    os.system(f'mkdir {working_dir}/ss')
except:
    pass
working_dir = f'{working_dir}/ss'

reference_dir = f'{data_dir}/reference/{args.population}'

simulation_size = args.size    # not currently used for anything 


#######################################

formatting_cmds = []
shapeit5_cmds = []
plink2_cmds = []
germline_cmds = []

# ALL IN ONE COMMAND: unzip the qced file, recode qced vcf to have 0/1, re-zip the 0/1 vcf, index the zipped vcf 

formatting_cmd = f"gunzip -c {vcf_file} > {working_dir}/{filehandle}.vcf; " + f"awk 'BEGIN {{OFS=\"\\t\"}} /^#/ {{print}} !/^#/ {{$4=\"0\"; $5=\"1\"; print}}' {working_dir}/{filehandle}.vcf > {working_dir}/{filehandle}_01.vcf; " + f"bcftools +fill-AN-AC {working_dir}/{filehandle}_01.vcf -Oz -o {working_dir}/{filehandle}_01.vcf.gz; " + f"bcftools index {working_dir}/{filehandle}_01.vcf.gz" 


formatting_cmds.append(formatting_cmd)

for chrom in range(1,23):

    # run shapeit5
    shapeit5_cmd = f"{shapeit5} -I {working_dir}/{filehandle}_01.vcf.gz -H {reference_dir}/chr{chrom}_01.vcf.gz -M {gmap_dir}/chr{chrom}.b38.gmap.gz -T 10 -O {working_dir}/{filehandle}_phased_chr{chrom}.bcf --log {working_dir}/{filehandle}_phased_chr{chrom}_shapeit.log --region {chrom}"
    shapeit5_cmds.append(shapeit5_cmd)

    # convert shapeit5 bcf to gzvcf -- incorporate this into the plink2 cmd actually 

    # run plink2 to generate .haps and .sample files
    plink2_cmd = f'bcftools view {working_dir}/{filehandle}_phased_chr{chrom}.bcf > {working_dir}/{filehandle}_phased_chr{chrom}.vcf.gz; ' +  f'{plink2} --vcf {working_dir}/{filehandle}_phased_chr{chrom}.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --export haps --out {working_dir}/{filehandle}_phased_chr{chrom}'
    plink2_cmds.append(plink2_cmd)

    # run germline
    germline_cmd = f'{germline2} {working_dir}/{filehandle}_phased_chr{chrom}.haps {working_dir}/{filehandle}_phased_chr{chrom}.sample {gmap_dir}/chr{chrom}_noheader.b38.gmap {working_dir}/{filehandle}_phased_chr{chrom}.match.hap -h'
    germline_cmds.append(germline_cmd)

    germline_cmd2 = f'{germline2} {working_dir}/{filehandle}_phased_chr{chrom}.haps {working_dir}/{filehandle}_phased_chr{chrom}.sample {gmap_dir}/chr{chrom}_noheader.b38.gmap {working_dir}/{filehandle}_phased_chr{chrom}.match'
    germline_cmds.append(germline_cmd2)

# execute steps in line 

print ('Running file formatting commands')
execute_commands(formatting_cmds, 20)

print ('Running SHAPEIT5 commands')
execute_commands(shapeit5_cmds, 20)

print ('Running haps/sample gen commands')
execute_commands(plink2_cmds, 20)

print ('Running germline2 commands')
execute_commands(germline_cmds, 20)

# Concatenate .match file results 
# look at every file in test_dir and if it ends in .match, write its entire contents to a new file 

output_file = f'{working_dir}/segment_input_temp.txt'
with open(output_file, 'w+') as f:
    pass


# Once for non-haploid output

match_files = [f.name for f in os.scandir(f'{working_dir}/') if f.name.endswith('.match')]
for f in match_files:
    chrom = f.split('chr')[-1].split('.')[0]
    with open(f'{working_dir}/{f}', 'r') as infile, open(f'{working_dir}/{f}2', 'w+') as outfile:
        for line in infile:
            line = line.strip('\n') + f'\t{chrom}\n' 
            outfile.write(line)

# concatenate
match2_files = [f.name for f in os.scandir(f'{working_dir}/') if f.name.endswith('.match2')]
for f in match2_files:
    with open(f'{working_dir}/{f}', 'r') as infile:
        with open(output_file, 'a') as outfile:
            outfile.write(infile.read())

# update .match file
match_df = pd.read_csv(output_file, header=None, sep='\t')
df_selected = match_df.iloc[:, [0,1,2,3,4,7]]

## FIX THE IDs BEING REPEATED 
df_selected[0] = df_selected[0].apply(lambda x: x.split('_')[-1])
df_selected[1] = df_selected[1].apply(lambda x: x.split('_')[-1])

df_selected.to_csv(f'{working_dir}/segment_input_final.txt', sep='\t', header=False, index=False)

print('Done with formatting')

try:
    os.system(f'mkdir {working_dir}/files')
except:
    pass

# move misc files out of top level

all_files = [x.name for x in os.scandir(working_dir) if x.is_file() and x.name != filehandlefull and x.name.split('_')[-1] != 'final.txt']
for f in all_files:
    os.system(f'mv {working_dir}/{f} {working_dir}/files/{f}')

print ('\nDone!\n')