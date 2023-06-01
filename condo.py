import re
import os
import argparse
import multiprocessing
import subprocess

def get_args():
  parser = argparse.ArgumentParser(description="Find the copy number of a feature in the genome")
  parser.add_argument("--no_uniq", action='store_true', default = False)
  #parser.add_argument("-n", "--name", help="The name of the feature(s)", default = "placeholder")
  parser.add_argument("-k", nargs='+', help="Kmer size(s). Default of 31", default = 31)
  parser.add_argument("-f", nargs='+', help="A consensus sequence or single occurence for your feature(s) of interest", required = True)
  parser.add_argument("-bed", nargs='+', help="The bed file(s) with coordinates for your feature in the reference genome. Required if no_uniq is not set.")
  parser.add_argument("-r", help="The directory containing reads to be CN assessed", required = True)
  parser.add_argument("-g", nargs='+', help="The reference genome to use", required = True)
  parser.add_argument("--gzip", action='store_true', default=False)
  return parser.parse_args()

args = get_args()
#name = args.n
k = args.k
feature_ref = args.f
reads = args.r
genome = args.g
if args.no_uniq == True:
    no_uniq = True
else:
    no_uniq = False
    # if args.bed is not None:
    #     fbed = args.bed
    # else:
    #     print("You must provide a bedfile with feature coordinates unless the no_uniq paramter is set! This enables finding unique feature kmers.")
    #     fail
fbed = args.bed
if args.gzip ==True:
    gzip = True
else:
    gzip = False

os.system('mkdir jellyfish_files results matched_windows feature')


#Get the sample IDs from the read names, propogate into config for rest of pipeline
command = 'ls ' + str(reads) + '| grep ".fq\|.fastq\|.fasta\|.fa"'
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)
#print(p.communicate())
files = []
while True:
  line = p.stdout.readline()
  if line != "":
      files.append(line)
  if not line:
    break

#pout = p.communicate()[0]
#print(pout.split("''"))

#all_files = os.popen('ls' + string(reads) + '| grep .fq\|.fastq\|.fasta\|.fa').readlines()
#all_files = all_files.strip()


#Identify ID names from file names
id_list = []
for item in files:
    line = item.strip()
    print(line)
    regex = re.search(r'([\S,_]+)_[0-9]\.(fa|fasta|fq|fastq|fastq.gz|fasta.gz|fa.gz|fq.gz)',line)
    id = regex.group(1)
    id_list.append(id)
#take only unique file names
id_list = set(id_list)
id_list = list(id_list)


#Do the same but for feature IDs
print(feature_ref)
print(args.f)
#if no_uniq == False:
feature_list = []
for item in feature_ref:
    line = item.strip()
    print(line)
    regex = re.search(r'/*([a-z,A-Z,0-9,-,.,_]+)\.(fa|fasta)',line)
    feature = regex.group(1)
    feature_list.append(feature)
#take only unique file names
feature_list = set(feature_list)
feature_list = list(feature_list)



# #Do the same but for feature IDs
# if no_uniq == False:
#     feature_list = []
#     for item in fbed:
#         line = item.strip()
#         print(line)
#         regex = re.search(r'/*([a-z,A-Z,0-9,-,.,_]+)\.(bed)',line)
#         feature = regex.group(1)
#         feature_list.append(feature)
#     #take only unique file names
#     feature_list = set(feature_list)
#     feature_list = list(feature_list)


#create a config file with contents adjusted by the argparse options
with open("config.yml","a+") as fil:
    fil.write("#config file for CONDO snakemake pipeline using user parameters")
    fil.write("\n")
    fil.write("FPATH:")
    fil.write("\n")
    for feature in feature_ref:
        fil.write("  - \"" + feature + "\"")
        fil.write("\n")
    fil.write("K:")
    fil.write("\n")
    for item in k:
        fil.write("  - \"" + str(item) + "\"")
        fil.write("\n")
    fil.write("ID:")
    fil.write("\n")
    for id in id_list:
        fil.write("  - \"" + id + "\"")
        fil.write("\n")
    fil.write("GENOME:")
    fil.write("\n")
    for path in genome:
        fil.write("  - \"" + path + "\"")
        fil.write("\n")
    fil.write("UNIQUE:")
    fil.write("\n")
    if no_uniq == False:
        fil.write("  - \"" + "unique" + "\"")
        fil.write("\n")
    if no_uniq == True:
        fil.write("  - \"" + "no_unique" + "\"")
        fil.write("\n")
    fil.write("BED:")
    fil.write("\n")
    for bed in fbed:
        fil.write("  - \"" + bed + "\"")
        fil.write("\n")
    #if no_uniq == False:
    fil.write("FEATURE:")
    fil.write("\n")
    for feature in feature_list:
        fil.write("  - \"" + feature + "\"")
        fil.write("\n")
    fil.write("SEQ_DIR:")
    fil.write("\n")
    fil.write("  - \"" + reads + "\"")
    fil.write("\n")
    fil.write("GZIP:")
    fil.write("\n")
    if gzip == True:
        fil.write("  - \"" + "True" + "\"")
    if gzip == False:
        fil.write("  - \"" + "False" + "\"")
    fil.write("\n")

for n in range(len(feature_list)):
    item = feature_list[n]
    path = feature_ref[n]
    bed_path = fbed[n]
    command = 'mkdir ' + 'feature/' + str(item)
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)
    command = 'scp '+ path +  ' feature/' + str(item) + "/"
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)
    command = 'scp '+ bed_path +  ' feature/' + str(item) + "/"
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)

os.system('snakemake all')
