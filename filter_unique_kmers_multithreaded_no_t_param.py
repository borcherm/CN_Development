import re
import argparse
import numpy
import random
import math
import multiprocessing as mp
import os

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-roi", help="fasta file for the region of interest in the genome")
  #parser.add_argument("-bed_comp", help= "A bed file with the genomic complement of the region of interest (inverse)")
  parser.add_argument("-k")
  parser.add_argument("-rest", help="A fasta file for the whole genome but excluding the region of interest")
  # parser.add_argument("-copycutoff")
  # parser.add_argument("-less")
  return parser.parse_args()

args = get_args()
roi = args.roi
#bed_comp = args.bed_comp
k = int(args.k)
rest = args.rest
# moreless = args.less
# cutoff = args.copycutoff

def start_process():
    print('Starting', mp.current_process().name)

def gc_content(sequence):
    AT = 0
    GC = 0
    for character in sequence:
        if character == "A":
            AT += 1
        elif character == "T":
            AT += 1
        elif character == "G":
            GC += 1
        elif character == "C":
            GC += 1
        else:
            print("Non-uppercase ATGC character detected")
            break
    total = GC+AT
    ratio = GC/total
    return(round(ratio,2))

def reverse_complement(sequence):
    new_seq = ""
    for character in sequence:
        if character == "A":
            new_seq += "T"
        elif character == "T":
            new_seq += "A"
        elif character == "G":
            new_seq += "C"
        elif character == "C":
            new_seq += "G"
        else:
            print("Non-uppercase ATGC character detected")
            break
    reversed_seq = new_seq[::-1]
    return(reversed_seq)


def find_unique(assembly):
    print("Identifying feature-specific kmers")
    with open(assembly, "r+") as ref:
        for line in ref:
            seq = line.strip()
            if seq[0] == ">":
                continue
            for n in range(len(seq)-k+1):
                kmer = seq[n:n+k]
                kmer_revcomp = reverse_complement(kmer)
                if kmer in kmer_dic:
                    #print(kmer)
                    del kmer_dic[kmer]
                elif kmer_revcomp in kmer_dic:
                    del kmer_dic[kmer_revcomp]
    return(kmer_dic)




print("Generating array kmers")
with open(roi, "r+") as fh:
    kmer_dic = {}
    while True:
        seqname = fh.readline().strip()
        if seqname == '':
            break
        seq = fh.readline().strip()
        #seqID = re.search(r'\>([\S]+)\:[\S]+',seqname)
        #chrID = seqID.group(1)
        Array_GC = gc_content(seq)
        for n in range(len(seq)-k+1):
            kmer = seq[n:n+k]
            if kmer in kmer_dic:
                kmer_dic[kmer] += 1
            else:
                kmer_dic[kmer] = 1


#combine any kmers from the array set that are reverse complements
print("The length before deleting duplicates:")
print(len(kmer_dic))
duplicate_rev_comp_list = []
duplicates = False
for kmer in kmer_dic:
    kmer_revcomp = reverse_complement(kmer)
    if kmer_revcomp in kmer_dic:
        duplicates = True
        if kmer not in duplicate_rev_comp_list:
            kmer_dic[kmer] += kmer_dic[kmer_revcomp]
            duplicate_rev_comp_list.append(kmer_revcomp)
if duplicates == True:
    print("There are reverse complement kmers present in the genomic sequence for this feature")
for entry in duplicate_rev_comp_list:
    del kmer_dic[entry]
print("The length after deleting duplicates:")
print(len(kmer_dic))



# below_cutoff = []
# if moreless == "more":
#     for item in kmer_dic:
#         if int(kmer_dic[item]) < int(cutoff):
#             below_cutoff.append(item)
# elif moreless == "less":
#     for item in kmer_dic:
#         if int(kmer_dic[item]) > int(cutoff):
#             below_cutoff.append(item)
# for kmer in below_cutoff:
#     del kmer_dic[kmer]



#split the assembly into single chromosomes, loop over each chromosome as a separate file
file_list = []
with open(rest, "r+") as ref:
    while True:
        seqname = ref.readline().strip()
        if seqname == "":
            break
        seq = ref.readline().strip()
        seqID = re.search(r'\>([\S]+)',seqname)
        if seqID.group(1) == False:
            print(seqname)
        chrmID = seqID.group(1)
        file_string = chrmID + "_assembly_split.fa"
        file_list.append(file_string)
        out = open(file_string,"w+")
        out.write(seqname)
        out.write("\n")
        out.write(seq)
        out.write("\n")
        out.close()


if __name__ == '__main__':
    pool = mp.Pool(processes= len(file_list),initializer=start_process)
    pool_outputs1 = pool.map(find_unique, [fh for fh in file_list])
    pool.close()
    pool.join()

print("The number of pool outputs is:")
print(len(pool_outputs1))
num_outputs = len(pool_outputs1)

kmer_dic_list=numpy.empty(num_outputs,dtype='object')
for i in range(num_outputs):
    kmer_dic_list[i] = pool_outputs1[i]


#merge the 23+ output dictionaries, removing keys which were not present in any single output
#these represent non unique kmers
for i in range(len(kmer_dic_list)):
    if i == 0:
        merged_kmer_dic = kmer_dic_list[i]
        print(len(merged_kmer_dic))
    else:
        #To allow myself to delete kmers during loop, I have to delete from a copy dictionary
        #This dicitonary becomes the merged dic at the end
        for item in merged_kmer_dic.copy():
            if item in kmer_dic_list[i].keys():
                continue
            elif kmer_revcomp in kmer_dic_list[i].keys():
                continue
            else:
                del merged_kmer_dic[item]

print(merged_kmer_dic)


# graph = open("kmer_frequency_x_array.tsv","a+")
# for entry in kmer_dic:
#     graph.write(entry+"\t"+str(kmer_dic[entry]))
#     graph.write("\n")

set = open("region_specific_kmers.fa","w+")
for entry in merged_kmer_dic:
    set.write(">" + str(merged_kmer_dic[entry]))
    set.write("\n")
    set.write(entry)
    set.write("\n")
set.close()

for file in file_list:
    os.remove(file)
