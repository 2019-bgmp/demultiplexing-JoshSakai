#!/usr/bin/env python

import numpy as np
import gzip
import matplotlib.pyplot as mpl

import argparse

def get_args():
    """inputting the file names for read 1, read 2, index 1, and index 2 """

    parser = argparse.ArgumentParser(description=' File name to be used')
    parser.add_argument('-r1', '--read1', help='Read 1 fastq file, must be zipped', type=str, required=True)
    parser.add_argument('-r2', '--read2', help='Read 2 fastq file, must be zipped', type=str, required=True)
    parser.add_argument('-i1', '--index1', help='Index 1 fastq file, must be zipped', type=str, required=True)
    parser.add_argument('-i2', '--index2', help='Index 2 fastq file, must be zipped', type=str, required=True)
    parser.add_argument('-q', '--quality', help='phred33 quality score threshold for reads', type=int, required=True)
    parser.add_argument('-b', '--barcodes', help='Known barcodes file as tsv with code and index on each line', type=str, required=True)
        #NOTE, barcodes file must conform to formating in the help documentation
    return parser.parse_args()

args= get_args()
r1 = args.read1
r2 = args.read2
i1 = args.index1
i2 = args.index2
q = args.quality
b = args.barcodes

def header_edit(header,i1,i2):
    """appends the given 'index1-index2' to the end of the header line and then returns a value that contains the edited header"""
    bars = ' '+i1+'-'+i2
    new_header = header + bars
    return new_header

def rv_comp(seq):
    """converts the DNA string into it's reverse compliment """
    rv = seq[::-1]
    flip = ''
    for initial in rv:
        if initial == 'A':
            flip += 'T'
        elif initial == 'T':
            flip += 'A'
        elif initial == 'C':
            flip += 'G'
        elif initial == 'G':
            flip += 'C'
        elif initial == 'N':
            flip += 'N'

    return flip

def low_Q(q, qual):
    """modified function that converts a single character into a phred score"""
    """the mean quality phred score is compared to the inputted quality threshold value and returns true if it is equal or greater and false if it is not"""
    value=0
    for cha in qual:
        #converts phred score
        score = (((ord(cha) - 33)))
        value += score
    if (int(value)/len(qual))>=q:
        #sums the quality scores and compares them to the desired value given during the argparse
        return True
    else:
        return False

def barcode_ref(b):
    """given the input file from the argparse this takes the barcodes within it and creates an array of barcodes as a key for demultiplexing"""
    """importantly the barcode file must be tab seperated and have a barcode on each line"""
    with open(b,'r') as b_ref :
        bar_list = []
        for line in b_ref:
            x = line.split('\t')[1]
            y = x.strip()
            bar_list.append(y)
        return bar_list
b_ref = barcode_ref(b)


def demult(r1,r2,i1,i2,q,b_ref):
    """Opening each filed passed in through the argparse with the zip function a set of temporary arrays and dictionaries are created to store each set of reads from the 4 files individually"""
    """this will increase the speed of this program since the sequencing files need to be only opened once while this program is running"""
    """Once a set of 4 lines have been read into the arrays then a reverse complement of the index2 sequence line is taken and then all quality scores for each read are compared to their given cutoff values.
    If all reads pass this benchmark the reverse complement is then compared to the index1 sequence line and if they are equivalent then the r1 and r2 inputs are then
    appended to a fw/rv file with the index1 prefix in the title."""
    """If the mean read quality for any of the reads does not match the set value then it is outputed into the lowq fw/rv outs"""
    """If the reads do meet the quality requirement but have indexes that do not match the b_ref sets then they are outputted to the index_hopped fw/rv outs"""
    """If a barcode does not match any of the b_ref codes, ie contains N or has a sequencing error, it will also be outputted into the lowq fw/rv outs"""

    with gzip.open(r1, 'rt') as r1, gzip.open(r2, 'rt') as r2, \
    gzip.open(i1, 'rt') as i1, gzip.open(i2, 'rt') as i2, \
    open('lowq_R1.fq', 'a') as bad_f, open('lowq_R2.fq', 'a') as bad_r, \
    open('index_hopped_R1.fq', 'a') as ind_hop_f, open('index_hopped_R2.fq', 'a') as ind_hop_r:
        codes = []
        for barcode in b_ref:
            codes.append(barcode + '_R1.fq')
            codes.append(barcode + '_R2.fq')
        files = {code: open(code,'a') for code in codes}
        counts = {'lowq':0, 'index_hopped':0}
        start = {}
        swap = {}
        #initializing a set of nested dictionaries to store the
        for bar in b_ref:
            start[bar]=0
        for bar in b_ref:
            swap[bar]={}
            for group in start:
                swap[bar][group]=0

        ln = 0
        r1_temp = []
        r2_temp = []
        i1_temp = []
        i2_temp = []
        i2_rv = []
        for r1_set , r2_set, i1_set, i2_set in zip(r1,r2,i1,i2):
            ln += 1
            r1_temp.append(r1_set.strip())
            r2_temp.append(r2_set.strip())
            i1_temp.append(i1_set.strip())
            i2_temp.append(i2_set.strip())
            if ln % 4 == 2:
                i2_rv = rv_comp(i2_set.strip())
            if ln % 4 == 0:
                if low_Q(q,r1_temp[3]) and low_Q(q,r2_temp[3]) \
                and low_Q(q,i1_temp[3]) and low_Q(q,i2_temp[3]) is True:
                    if str(i1_temp[1]) == str(i2_rv) \
                    and str(i1_temp[1]) in b_ref: #if everything comes up good
                        header_r1 = header_edit(r1_temp[0], i1_temp[1], i2_temp[1])
                        header_r2 = header_edit(r2_temp[0], i1_temp[1], i2_temp[1])
                        files[(i1_temp[1] + "_R1.fq")].write(header_r1 + '\n' + r1_temp[1] + '\n' + r1_temp[2] + '\n' + r1_temp[3] + '\n')
                        files[(i1_temp[1] + "_R2.fq")].write(header_r2 + '\n' + r2_temp[1] + '\n' + r2_temp[2] + '\n' + r2_temp[3] + '\n')
                        if i1_temp[1] in counts: #if this barcode has been read before
                            counts[i1_temp[1]]+=1
                            r1_temp.clear()# clearing out the arrays so that the next read can be inputted
                            r2_temp.clear()
                            i1_temp.clear()
                            i2_temp.clear()
                        else: #if this barcode has not been read before
                            counts[i1_temp[1]] = 1
                            r1_temp.clear()# clearing out the arrays so that the next read can be inputted
                            r2_temp.clear()
                            i1_temp.clear()
                            i2_temp.clear()
                    elif str(i1_temp[1]) != str(i2_rv) \
                    and str(i1_temp[1]) in b_ref \
                    and str(i2_rv) in b_ref:# if the reverse complement of the index2 does not match index1 /if its index hopped
                        header_r1 = header_edit(r1_temp[0], i1_temp[1], i2_temp[1])
                        header_r2 = header_edit(r2_temp[0], i1_temp[1], i2_temp[1])
                        ind_hop_f.write(header_r1 + '\n' + r1_temp[1] + '\n' + r1_temp[2] + '\n' + r1_temp[3] + '\n')
                        ind_hop_r.write(header_r2 + '\n' + r2_temp[1] + '\n' + r2_temp[2] + '\n' + r2_temp[3] + '\n')
                        counts['index_hopped'] += 1
                        swap[i1_temp[1]][i2_rv] += 1
                        r1_temp.clear()# clearing out the arrays so that the next read can be inputted
                        r2_temp.clear()
                        i1_temp.clear()
                        i2_temp.clear()
                    else:#if it has N's in the barcode or is a sequencing error
                        header_r1 = header_edit(r1_temp[0], i1_temp[1], i2_temp[1])
                        header_r2 = header_edit(r2_temp[0], i1_temp[1], i2_temp[1])
                        bad_f.write(header_r1 + '\n' + r1_temp[1] + '\n' + r1_temp[2] + '\n' + r1_temp[3] + '\n')
                        bad_r.write(header_r2 + '\n' + r2_temp[1] + '\n' + r2_temp[2] + '\n' + r2_temp[3] + '\n')
                        counts['lowq'] += 1
                        r1_temp.clear()# clearing out the arrays so that the next read can be inputted
                        r2_temp.clear()
                        i1_temp.clear()
                        i2_temp.clear()
                else:#if the quality of the sequences is too low.
                    header_r1 = header_edit(r1_temp[0],i1_temp[1],i2_temp[1])
                    header_r2 = header_edit(r2_temp[0],i1_temp[1],i2_temp[1])
                    bad_f.write(header_r1 + '\n' + r1_temp[1] + '\n' + r1_temp[2] + '\n' + r1_temp[3] + '\n')
                    bad_r.write(header_r2 + '\n' + r2_temp[1] + '\n' + r2_temp[2] + '\n' + r2_temp[3] + '\n')
                    counts['lowq'] += 1
                    r1_temp.clear()# clearing out the arrays so that the next read can be inputted
                    r2_temp.clear()
                    i1_temp.clear()
                    i2_temp.clear()
    for file in files.values():
        file.close()
    swapped = swap
    """taking our group values this creates the values that will go into our outfile containing the number of reads and % of total reads each category comprises"""
    for group in counts:
        if group !='lowq' and group !='index_hopped':
            swap[group][group] = counts[group]
    print ("Barcode\tCount\t% of Reads")
    for group in counts:
        print (str(group)+'\t'+str(counts[group]) + '\t' + str(counts[group]/(ln/4)))
    return counts, swap, swapped

"""a set of graphing functions to display the demultiplexing run results"""

def countplot(counts):
    y = []
    x = []
    for group in counts:
        x.append(group)
        y.append(counts[group])
    mpl.bar(x, y, align='center', alpha=0.5)
    mpl.ylabel('Count')
    mpl.xlabel('Barcode')
    mpl.xticks(np.arange(26), x, rotation=90)
    mpl.title('Counts for each barcode')
    mpl.savefig('counts.png', bbox_inches = "tight")
    mpl.close()
    return

def swapplot(swap):
    z = []
    twod = []
    for d in swap:
        for k in swap[d]:
            if len(z) < 24:
                z.append(swap[d][k])
            else:
                twod.append(z)
                z = [swap[d][k]]
    b_ref = barcode_ref(b)
    x = range(24)
    y = range(24)
    intensity = np.array(twod)
    x,y = np.meshgrid(x,y)
    mpl.pcolormesh(x, y, intensity)
    mpl.colorbar()
    mpl.xticks(np.arange(24),b_ref,rotation=90)
    mpl.yticks(np.arange(24),b_ref)
    mpl.title('Counts for each combination of barcodes')
    mpl.gcf().subplots_adjust(left=0.2)
    mpl.savefig('all_count.png', bbox_inches="tight")
    mpl.close()
    return

def swapedplot(swapped):
    z = []
    twod = []
    for d in swapped:
        for k in swapped[d]:
            if len(z) < 24:
                z.append(swapped[d][k])
            else:
                twod.append(z)
                z = [swapped[d][k]]
    b_ref = barcode_ref(b)
    x = range(24)
    y = range(24)
    intensity = np.array(twod)
    x,y = np.meshgrid(x,y)
    mpl.pcolormesh(x, y, intensity)
    mpl.colorbar()
    mpl.xticks(np.arange(24),b_ref,rotation=90)
    mpl.yticks(np.arange(24),b_ref)
    mpl.gcf().subplots_adjust(left=0.2)
    mpl.title('Counts for each combination of barcodes excluding paired barcodes')
    mpl.savefig('no_pair_count.png',bbox_inches="tight")
    mpl.close()
    return

counts, swap, swapped = demult(r1, r2, i1, i2, q, b_ref)

countplot(counts)
swapplot(swap)
swapedplot(swapped)

"""outfiles were then gziped from the command line afterwards for ease of submission"""
