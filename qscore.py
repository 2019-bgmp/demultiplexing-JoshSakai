#!/usr/bin/env python3

import argparse
import numpy as np
import gzip

def start_args():
    parser=argparse.ArgumentParser("Script to create a distribution of quality scores of fastq paired end reads")
    parser.add_argument("-f", "--input_file", help="Input file", required=True)
    parser.add_argument("-l", "--length", help="Target line length, needed to construct list parameters", type=int, required=True)
    parser.add_argument("-o", "--output_file", help="output file", required=True)
    return parser.parse_args()
    
args = start_args()             #grabs the flags for the global statement
file = args.input_file                #sets the name of the input file in my naming convention of choice
output = args.output_file                #sets the name of the output file in my naming convention of choice
target = args.length


def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter)-33

all_qscores = np.zeros((target),dtype=float)

with gzip.open(file,"rt") as base:
    LN=0
    num_records=0
    for line in base:
        LN+=1
        if LN%4 == 0:
            num_records+=1
            line = line.strip('\n')
            cha=0
            for letter in line:
                q_score = convert_phred(letter)
                all_qscores[cha]+=q_score
                cha+=1
print(all_qscores)
print(num_records)

position = []
index_value=0
print("# Base Pair", "\t", "Mean Quality Score")
for avg in all_qscores:
    all_qscores[index_value]=avg/num_records
    print(index_value, "\t", all_qscores[index_value])
    position.append(index_value)
    index_value += 1
    

import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
plt.xlabel('Base pair value')
plt.ylabel('mean quality score')
plt.title('Mean Quality Score per Base')
plt.bar(position, all_qscores)
plt.savefig(output+'_graph.png')