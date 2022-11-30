import re
import argparse
import numpy
import random
import math
import statistics
import multiprocessing as mp
import gzip
from collections import Counter

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-r1")
  parser.add_argument("-r2", required = False)
  parser.add_argument("-nc1")
  parser.add_argument("-nc2", required = False)
  parser.add_argument("-ngcn")
  parser.add_argument("-ID")
  parser.add_argument("-gcn")
  return parser.parse_args()

def takeSecond(elem):
    return int(elem[1])


args = get_args()
counts_r1 = args.r1
if args.r2 is not None:
    counts_r2 = args.r2
else:
    counts_r2 = False
norm_counts1 = args.nc1
if args.nc2 is not None:
    norm_counts2 = args.nc2
else:
    norm_counts2 = False
norm_gcn = args.ngcn
ID = args.ID
gcn= args.gcn


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
    revcomp = reverse_complement(kmer)
    if kmer in genomic_copy_number_dic:
        genomic_cn_normalized_counts_dic[kmer] = kmer_dic[kmer]/genomic_copy_number_dic[kmer]
    elif revcomp in genomic_copy_number_dic:
        genomic_cn_normalized_counts_dic[kmer] = kmer_dic[kmer]/genomic_copy_number_dic[revcomp]
    else:
        fail
    genomic_cn_normalized_counts_list.append(genomic_cn_normalized_counts_dic[kmer])

print(genomic_cn_normalized_counts_dic)

#find median
median = int(statistics.median_low(genomic_cn_normalized_counts_list))
print(median)

mean = int(statistics.mean(genomic_cn_normalized_counts_list))
print(mean)

norm_counts_dic = Counter()
with open(norm_counts1, "r+") as norm1:
    while True:
        l1 = norm1.readline().strip()
        kmer = norm1.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        norm_counts_dic[kmer] += count

if norm_counts2 is not False:
    with open(norm_counts2, "r+") as norm2:
        while True:
            l1 = norm2.readline().strip()
            kmer = norm2.readline().strip()
            if l1 == "":
                break
            regex = re.search(r'\>([\S]+)',l1)
            count = int(regex.group(1))
            norm_counts_dic[kmer] += count

norm_gcn_dic = Counter()
with open(norm_gcn) as ngcn:
    while True:
        l1 = ngcn.readline().strip()
        kmer = ngcn.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        norm_gcn_dic[kmer] += count

norm_counts_adj = Counter()
for item in norm_counts_dic:
    revcomp = reverse_complement(item)
    if item in norm_gcn_dic:
        norm_counts_adj[item] = float(norm_counts_dic[item]/norm_gcn_dic[item])
    elif revcomp in norm_gcn_dic:
        norm_counts_adj[item] = float(norm_counts_dic[item]/norm_gcn_dic[revcomp])
    else:
        fail


norm_median = statistics.median(list(norm_counts_adj.values()))
print(norm_median)
diploid_norm_median = norm_median/2

ribosome_copy_number = median/norm_median
print("Median Copy Number" + "\n")
print(ribosome_copy_number)

mean_copy_number = mean/norm_median
print("Mean Copy Number:" + "\n")
print(mean_copy_number)

diploid_ribosome_copy_number = median/diploid_norm_median
print("Diploid Median Copy Number" + "\n")
print(diploid_ribosome_copy_number)

diploid_mean_copy_number = mean/diploid_norm_median
print("Diploid Mean Copy Number:" + "\n")
print(diploid_mean_copy_number)

with open("Copy_Numbers.tsv", "a+") as summary:
    #summary.write("Sample" + "\t" + "Median_Haploid_Copy_Number" + "\t" + "Median_Diploid_Copy_Numer" + "\t" + "Mean_Haploid_Copy_Number" + "\t" + "Mean_Diploid_Copy_Number" + "\n")
    summary.write(ID + "\t" +str(round(ribosome_copy_number)) + "\t" + str(round(diploid_ribosome_copy_number)) + "\t" + str(round(mean_copy_number)) + "\t" + str(round(diploid_mean_copy_number))+"\n")
