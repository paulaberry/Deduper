#!/usr/bin/env python3
import gzip
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import os
import re

# We need to remove PCR duplicates. The duplicates will have these characteristics:
#   1. Same chromosome
#     --found with RNAME (SAM col 3)
#   2. Same strand
#       --found with FLAG (SAM col 2)
#   3. Same starting position
#     --found with POS (SAM col 4)
#     --Soft clipping can make this hard to figure out!
#   4. Same UMI
#     --found with QNAME (SAM col 1)
#
# Pitfalls:
#   Samtools sorts by QNAME but not by UMI. So we will sort by chromosomal coordinates which means all entries that occur on the same chromosome will have to be checked against all others that occur on the same chromosome.

def get_args():
    """Function to pass in the SAM file and options to remove PCR duplicates."""
    getfiles = argparse.ArgumentParser(description="A program to remove PCR duplicates from SAM files.")
    getfiles.add_argument("-f", "--samfile", help = "SAM file to be deduped", required = True)
    getfiles.add_argument("-p", "--paired", action = "store_true", help = "reads are paired-end if option is selected", required = False)
    getfiles.add_argument("-u", "--UMIfile", help = "text file containing a list of UMIs", required = False, default = "None")
    getfiles.add_argument("-q", "--quality", action = "store_true", help = "if present, keep duplicate read with higher quality score", required = False)
    getfiles.add_argument("-w", "--write", action = "store_true", help = "if present, duplicates and discarded reads are written to a separate file", required = False)
    getfiles.add_argument("-r", "--randomers", help = "if present, UMIs are randomers, incompatible with --UMIfile option", required = False, default = 0)
    return getfiles.parse_args()
args = get_args()

# required arguments
filename = str(args.samfile)

# optional arguments with strings
umi_file = str(args.UMIfile)

# optional boolean arguments
pairedend = bool(args.paired)
qualitysort = bool(args.quality)
writeout = bool(args.write)

# optional integer arguements
umi_length = int(args.randomers)

# Warnings for incompatible arguments
if umi_file != "None" and umi_length != 0:
    print("Choose randomer length or specify UMI list file!")
elif umi_file == "None" and umi_length == 0:
    print("Choose randomer length or specify UMI list file!")

# functions
def strand(flag):
    """A function that takes a bitwise FLAG (column 2) as an argument, and returns the strand."""
    strand = "+" #set default direction to fwd
    if ((int(flag) & 16) == 16): # eval bitwise flag
        strand = "-"
    return strand

def findposition(read_direction, pos, CIGAR):
    """A function that takes as arguments the direction, position, and CIGAR strings from a SAM file entry.

    If "S" is the first letter (not character!) in the string we will have to adjust the read for soft-clipping when read direction is forward. If "S" is the last character in the string, this function will adjust the read position for soft clipping based on the last set of numbers when the read direction is reverse. If the read direction is reverse, the CIGAR string will be tallied and added to the number in column 3 of the SAM file, and if the resulting numbers are the same, the positions are the same.
        Input from SAM col 2, 3, 5: (0, 15, 2S15M, 13, 17M, None, None, None, None) returns True.
        Input from SAM col 2, 3, 5: (0, 15, 17M, 15, 2S15M,  None, None, None, None) returns False.
        Input from SAM col 2, 3, 5: (16, 150, 8M3S, 130, 20S3M1I3M1D5M, 15, 13M, 15, 13M) returns True.
        Input from SAM col 2, 3, 5: (16, 150, 17M2S, 150, 19M, 40, 20M, 40, 20M) returns False."""
    S = 0 # Softclip at beginning of fwd reads is subtracted, at end of reverse reads is added to position
    D = 0 # Deleted nts in reverse reads is added to position
    # Inserted nts in reverse reads is null to position, do not need
    N = 0 # Skipped regions in reverse reads is added to position
    M = 0 # Matched nts in reverse reads is added to position
    if read_direction == "+":
        softclip = re.findall(r"^\d+S", CIGAR)
        if softclip != []:
            S = int(softclip[0][:-1])
        if softclip != []:
            softclip = int(softclip[0][:-1])
        pos = int(pos) - S
        return pos

    elif read_direction == "-":
        softclip = re.findall(r"\d+S$", CIGAR)
        delete = re.findall(r"\d+D", CIGAR)
        skip = re.findall(r"\d+N", CIGAR)
        match = re.findall(r"\d+M", CIGAR)
        if softclip != []:
            S = int(softclip[0][:-1])
        if delete != []:
            delcount = 0
            for i in delete:
                D = D + int(delete[delcount][:-1])
                delcount = delcount + 1
        if skip != []:
            skipcount = 0
            for i in skip:
                N = N + int(skip[skipcount][:-1])
                skipcount = skipcount + 1
        if match != []:
            matchcount = 0
            for i in match:
                M = M + int(match[matchcount][:-1])
                matchcount = matchcount + 1
        pos = int(pos) + S + D + N + M
        return pos

# algorithm
# initialize UMI dictionary if -u is selected
if umi_file != "None":
    UMI_dict = {}
    count  = 0
    umi_list = open(umi_file, "r")
    for line in umi_list:
        count = count + 1
        line = line.strip()
        umi_length = len(line) # set variable for length of UMI to later extract from QNAME field
        UMI_dict[line] = count
    umi_list.close()

# Initialize files for writing out duplicated and unknown reads
if writeout == True:
    dupe_filename = "duplicate." + str(filename)
    dupe_file = open(dupe_filename, "w")
    dupe_file.seek(0) # make sure to start at beginning of file
    unknown_filename = "unknown." + str(filename)
    unknown_file = open(unknown_filename, "w")
    unknown_file.seek(0) # make sure to start at the beginning of the file

# Initialize deduplicated read file
deduped_filename = "deduped." + str(filename)
deduped_file = open(deduped_filename, "w")
deduped_file.seek(0) # make sure starting at beginning of file

samfile = open(filename, "r") # make sure to start at beginning of SAM file

# initialize dictionary to keep track of duplicate reads and chosen reads
deduped_dict = {}
duptup = "initialize" # initialize duptup

if pairedend == False: # algorithm for non-paired end reads
    while True:
        samline = samfile.readline().strip()
        if samline == "": # stop while loop at end of SAM file
            #print(deduped_dict)
            for key in deduped_dict:
                print(deduped_dict[key][1], file = deduped_file) # print out everything in the deduplicate dictionary
            break
        elif samline.startswith("@"):
            print(str(samline), file = deduped_file) # print header lines to deduped file
        else:
            #print(samline)
            samlist = samline.split()
            UMI = samlist[0][-umi_length:]
            chromosome = samlist[2]
            read_direction = "+" # default read direction is + strand
            if ((int(samlist[1]) & 16) == 16):
                read_direction = "-" # set read direction if L1 is mapped to - strand
            position = findposition(read_direction, int(samlist[3]), samlist[5])

            if duptup == "initialize": # first SAM read option
                duptup = (chromosome, UMI, read_direction, position)
                deduped_dict[duptup] = (samlist[4], str(samline))
                print("first sam entry into dictionary")
            elif chromosome != duptup[0]: # if we change chromosomes
                for key in deduped_dict:
                    print(deduped_dict[key][1], file = deduped_file) # print out everything in the dedupolicate dictionary
                deduped_dict = {} # reset dictionary to empty
                duptup = (chromosome, UMI, read_direction, position)
                deduped_dict[duptup] = (samlist[4], str(samline)) # first entry of new chromosome into dictionary
            else:
                duptup = (chromosome, UMI, read_direction, position)
                print("else branch started")
                if duptup in deduped_dict:
                        if qualitysort != True:
                            print(samline, file = dupe_file)
                        else:
                            if samlist[4] > deduped_dict[duptup][0]:
                                if writeout == True:
                                    print(deduped_dict[duptup][1], file = dupe_file)
                                    deduped_dict[duptup] = (samlist[4], str(samline))
                                else:
                                    deduped_dict[duptup] = (samlist[4], str(samline))
                            else:
                                pass
                else:
                    if umi_file != "None":
                        if UMI not in UMI_dict:
                            if writeout == True:
                                print(samline, file= unknown_file)
                            else:
                                pass
                        else:
                            duptup = (chromosome, UMI, read_direction, position)
                            deduped_dict[duptup] = (samlist[4], str(samline))
                    else:
                        duptup = (chromosome, UMI, read_direction, position)
                        deduped_dict[duptup] = (samlist[4], str(samline))


elif pairedend == True:
    duptup = (-1, -1, -1, -1, -1) # initialize duptup
    while True:
        samlineR1 = samfile.readline().strip()
        samlineR2 = samfile.readline().strip()
        if samlineR1 == "": # stop while loop at end of SAM file
            break
        elif samlineR1.startswith("@"):
            print(str(samlineR1), file = deduped_file) # print header lines to deduped file
            if samlineR2.startswith("@"):
                print(str(samlineR2), file = deduped_file) # print header lines to deduped file
            else:
                samlineR1 = samlineR2
                samlineR2 = samfile.readline().strip()
                continue
        else:
            samlistR1 = samlineR1.split()
            samlistR2 = samlineR2.split()
            UMI = samlistR1[0][-umi_length:]
            chromosome = samlistR1[2]
            if chromosome != duptup[0]: # if we change chromosomes
                for i in deduped_dict:
                    print(deduped_dict[i][1], file = deduped_file)
                    print(deduped_dict[i][2], file = deduped_file) # print out everything in the dedupolicate dictionary
                deduped_dict = {} # reset dictionary to empty
            read_direction = "+" # default read direction is + strand
            if ((int(samlistR1[1]) & 16) == 16):
                read_direction = "-" # set read direction if L1 is mapped to - strand
            positionR1 = findposition(read_direction, int(samlistR1[3]), samlistR1[5])
            positionR2 = findposition(read_direction, int(samlistR2[3]), samlistR2[5])
            duptup = (chromosome, UMI, strand, positionR1, positionR2)
            if umi_file != "None" and writeout == True and UMI not in UMI_dict:
                print(samlineR1, file = unknown_file)
                print(samlineR2, file = unknown_file)
            if duptup in deduped_dict:
                if qualitysort != True:
                    print(samlineR1, file = dupe_file)
                    print(samlineR2, file = dupe_file)
                else:
                    if mean(samlistR1[10], samlistR2[10]) > deduped_dict[duptup][0]:
                        if writeout == True:
                            print(deduped_dict[duptup][1], file = dupe_file)
                        deduped_dict[duptup] = (samlistR1[10], str(samlineR1), samlistR2[10], str(samlineR2))
            else:
                deduped_dict[duptup] = (mean(samlistR1[4] , samlistR2[4]), str(samlineR1), str(samlineR2))
# close all the files
if writeout == True:
    dupe_file.close()
    unknown_file.close()
deduped_file.close()
samfile.close()
