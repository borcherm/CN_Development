import re
import argparse
import numpy
import random
import math
import statistics
import multiprocessing as mp
import gzip
import os
from collections import Counter

def get_args():
  parser = argparse.ArgumentParser(description="Find copy number of a feature in WGS data with pre-counted kmer multiplicities in the feature and a normalization set")
  parser.add_argument("-r1", help="A jellyfish counts fasta of kmers from the feature in the first WGS read (for paired and unpaired)")
  parser.add_argument("-r2", required = False, help="A jellyfish counts fasta of kmers from the feature in the second WGS read (if paired)")
  parser.add_argument("-ID", help="A sample ID to use in the output")
  parser.add_argument("-gcn", help="A jellyfish counts file of kmers from the feature of interest, so counts are indicating the multiplicity of the kmers in the feature.")
  parser.add_argument("-mode1", help="The mode kmer count of all kmers in the sequencing data read 1")
  parser.add_argument("-mode2")
  return parser.parse_args()

def takeSecond(elem):
    return int(elem[1])


args = get_args()
counts_r1 = args.r1
if args.r2 is not None:
    counts_r2 = args.r2
else:
    counts_r2 = False
ID = args.ID
gcn= args.gcn
read_mode_1 = int(args.mode1)
read_mode_2 = int(args.mode2)

counts_list = []
kmer_list = []
kmer_dic = {}

with open(counts_r1, "r+") as set:
    while True:
        l1 = set.readline().strip()
        kmer = set.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        kmer_dic[kmer] = count
        #kmer_list.append([kmer,count])

if counts_r2 is not False:
    with open(counts_r2, "r+") as set:
        while True:
            l1 = set.readline().strip()
            kmer = set.readline().strip()
            if l1 == "":
                break
            regex = re.search(r'\>([\S]+)',l1)
            count = int(regex.group(1))
            if kmer in kmer_dic:
                kmer_dic[kmer] = kmer_dic[kmer] + count
            else:
                kmer_dic[kmer] = count
#kmer_list.sort(key=takeSecond, reverse=True)
#print(kmer_list[0:100])
#print(kmer_list[-100:])

with open(gcn, "r+") as fh:
    genomic_copy_number_dic = {}
    while True:
        l1 = fh.readline().strip()
        kmer = fh.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        genomic_copy_number_dic[kmer] = count


genomic_cn_normalized_counts_dic = {}
genomic_cn_normalized_counts_list = []
for kmer in kmer_dic:
    genomic_cn_normalized_counts_dic[kmer] = kmer_dic[kmer]/genomic_copy_number_dic[kmer]
    genomic_cn_normalized_counts_list.append(genomic_cn_normalized_counts_dic[kmer])

print(genomic_cn_normalized_counts_dic)


standard_dev = statistics.stdev(genomic_cn_normalized_counts_list)
mean = statistics.mean(genomic_cn_normalized_counts_list)

final_list = []
for item in genomic_cn_normalized_counts_list:
    count = item
    if count != 0:
        if abs(mean-count) < 3*standard_dev:
            final_list.append(item)
        else:
            print(item)


#find median
median = int(statistics.median_low(final_list))
print("The median genomic cn adjusted feature kmer counts")
print(median)
mean = int(statistics.mean(final_list))
print("The mean genomic cn adjusted feature kmer counts")
print(mean)

genomic_norm_mode = round((read_mode_1 + read_mode_2)/2)*2
diploid_genomic_norm_mode = genomic_norm_mode/2
print(read_mode_1)
print(read_mode_2)
print(genomic_norm_mode)

ribosome_copy_number = median/genomic_norm_mode
print("Median Copy Number" + "\n")
print(ribosome_copy_number)

mean_copy_number = mean/genomic_norm_mode
print("Mean Copy Number:" + "\n")
print(mean_copy_number)

diploid_ribosome_copy_number = median/diploid_genomic_norm_mode
print("Diploid Median Copy Number" + "\n")
print(diploid_ribosome_copy_number)

diploid_mean_copy_number = mean/diploid_genomic_norm_mode
print("Diploid Mean Copy Number:" + "\n")
print(diploid_mean_copy_number)

exists = os.path.isfile("./Copy_Numbers.tsv")

summary = open("Copy_Numbers.tsv","a+")

if exists == False:
    summary.write("Sample" + "\t" + "Median_Haploid_Copy_Number" + "\t" + "Median_Diploid_Copy_Numer" + "\t" + "Mean_Haploid_Copy_Number" + "\t" + "Mean_Diploid_Copy_Number" + "\n")
    summary.write(ID + "\t" +str(round(ribosome_copy_number)) + "\t" + str(round(diploid_ribosome_copy_number)) + "\t" + str(round(mean_copy_number)) + "\t" + str(round(diploid_mean_copy_number))+"\n")
else:
    summary.write(ID + "\t" +str(round(ribosome_copy_number)) + "\t" + str(round(diploid_ribosome_copy_number)) + "\t" + str(round(mean_copy_number)) + "\t" + str(round(diploid_mean_copy_number))+"\n")
summary.close()
