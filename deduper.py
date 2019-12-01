#!/usr/bin/env python3
import gzip
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import os

def get_args():
    """Function to pass in the SAM file and options to remove PCR duplicates."""
    getfiles = argparse.ArgumentParser(description="A program to remove PCR duplicates from SAM files.")
    getfiles.add_argument("-f", "--samfile", help = "SAM file to be deduped", required = True)
    getfiles.add_argument("-p", "--paired", action = "store_true", help = "reads are paired-end if option is selected", required = False)
    getfiles.add_argument("-u", "--UMIfile", help = "text file containing a list of UMIs", required = False)
    getfiles.add_argument("-h", "--help", help = "A program to remove PCR duplicates from SAM files.", required = False)
    getfiles.add_argument("-q", "--quality", action = "store_true", help = "if present, keep duplicate read with higher quality score", required = False)
    getfiles.add_argument("-w", "--write", action = "store_true", help = "if present, duplicates and discarded reads are written to a separate file", required = False)
    getfile.add_argument("-d", "--destination", help = "if present, specifies output filepath", required = False)
    return getfiles.parse_args()
args = get_args()

# required arguments
filename = str(args.samfile)

# optional arguments with strings
output_path = str(args.destination)
umi_file = str(args.UMIfile)

# optional boolean arguments
pairedend = bool(args.paired)
qualitysort = bool(args.quality)
writeout = bool(args.write)

# functions
# Strategy: Create series of functions that return True or False, if all are True, discard the duplicate according to options set when intializing the program
def match_UMI(UMIL1R1, UMIL2R1, UMIL1R2 = None, UMIL2R2 = None):
    """A function that takes two UMIs (QNAME field from SAM file) as arguments and returns True if the two UMI's match, and False if they do not. R2 UMIs are set to None as default, and will only be used if a paired-end read SAM file is being processed.
        Input from SAM col 1: (AAAA, AAAA, None, None) returns True.
        Input from SAM col 1: (AAAA, AATT, None, None) returns False.
        Input from SAM col 1: (AAAA, AAAA, TTTT, TTTT) returns True.
        Input from SAM col 1: (AAAA, AATT, TTTT, TTAA) returns False."""
    if UMIL1R2 == None:
        if UMIL1R1 == UMIL2R1:
            return True
        else:
            return False
    else:
        if UMIL1R1 == UMIL2R1:
            if UMIL1R2 == UMIL2R2:
                return True
            else: return False
        else: return False

def match_strand(flagL1R1, flagL2R1, flagL1R2 = None, flagL2R2 = None):
    """A function that takes two bitwise FLAGs (column 2) as arguments, and returns True if they match and False if they do not.
        Input from SAM col 2: (16, 16, None, None) returns True.
        Input from SAM col 2: (16, 0, None, None) returns False.
        Input from SAM col 2: (16, 16, 16, 16) returns True.
        Input from SAM col 2: (16, 0, 16, 16) returns False."""
    if flagL1R2 == None:
        if flagL1R1 == flagL2R1:
            return True
        else:
            return False
    else:
        if flag L1R1 == flagL2R1:
            if flag L1R2 == flagL2R2:
                return True
            else:
                return False
        else:
            return False

def match_chromosome(chromosomeL1R1, chromosomeL2R1, chromosomeL1R2 = None, chromosomeL2R2 = None):
    """A function that takes two chromosome names (RNAME field from SAM file) as arguments and returns True if the two chromosome names match and False if they do not.
        Input from SAM col 3: (2, 2, None, None) returns True.
        Input from SAM col 3: (1, 2, None, None) returns False.
        Input from SAM col 3: (1, 1, 1, 1) returns True.
        Input from SAM col 3: (1, 2, 1, 2) returns False."""
    if chromosomeL1R2 == None:
        if chromosomeL1R1 == chromosomeL2R1:
            return True
        else:
            return False
    else:
        if chromosomeL1R1 == chromosomeL2R1:
            if chromosomeL1R2 == chromosomeL2R2:
                return True
            else:
                return False
        return False

def match_position(read_direction, positionL1R1, CIGARL1R1, positionL2R1, CIGARL2R1, positionL1R2 = None, CIGARL1R2 = None, positionL2R2 = None, CIGARL2R2 = None):
    """A function that takes as arguments the position and CIGAR strings from two lines in a SAM file. The function will account for softclipping from the CIGAR strings and update the values of the position variables. The function will then return True if the position variable values are the same, and False if they are not. For forward reads, only soft clipping at the beginning of the CIGAR string needs to be accounted for. For reverse reads, the entire CIGAR string must be tallied to find the actual position.

    If "S" is the first letter (not character!) in the string we will have to adjust the read for soft-clipping when read direction is forward. If "S" is the last character in the string, this function will adjust the read position for soft clipping based on the last set of numbers when the read direction is reverse. If the read direction is reverse, the CIGAR string will be tallied and added to the number in column 3 of the SAM file, and if the resulting numbers are the same, the positions are the same.
        Input from SAM col 2, 3, 5: (0, 15, 2S15M, 13, 17M, None, None, None, None) returns True.
        Input from SAM col 2, 3, 5: (0, 15, 17M, 15, 2S15M,  None, None, None, None) returns False.
        Input from SAM col 2, 3, 5: (16, 150, 8M3S, 130, 20S3M1I3M1D5M, 15, 13M, 15, 13M) returns True.
        Input from SAM col 2, 3, 5: (16, 150, 17M2S, 150, 19M, 40, 20M, 40, 20M) returns False."""

    pass
    #if positionL1R2 == None:



# algorithm

# initialize UMI dictionary if -u is selected
UMI_dict = {}
if umi_file in locals():
    f = open("umi_file", "r")
    count = 0
    for line in f:
        count = count + 1
        line = line.strip()
        UMI_dict[line] = count



# # problem
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

# FUNCTIONS

#
#

# Initialize samfile_deduped.sam
#
#
# # Deduping
# If -p is not present we will be going a line at a time using readline(). If -p is present, we will be going four lines at a time. The main algorithm will be written assuming -p is not specified.
#
# -u is present, and the QNAME does not match any keys in UMI_dict, it is written to unknown.sam and the next line read in.
#
# L1R1 = first line in set. This variable stays static until the line the program is comparing it to has a different Chromosome name (RNAME). It then becomes the line below, and the algorithm restarts from the top.
#
# L2R1 = second line in set, being compared to L1R1. This line changes as the program loops over all lines that have the same Chromosome name (RNAME).
#
#   if -p is specified:
# L1R2 = third line in set, same rules as for L1R1
# L2R2 = fourth line in set, same rules as for L2R1
#
# if line starts with "@" ignore, write to output file and move onto next line until both (or all four in the case of Paired-End reads) are valid SAM entries.
#
#
#   Split the strings into lists. List position will be as follows:
#     0 =	QNAME	String	Query template NAME
#     1 = FLAG	Int	bitwise FLAG
#     2 = RNAME	String	References sequence NAME
#     3 = POS	Int	1- based leftmost mapping POSition
#     4 = MAPQ	Int	MAPping Quality
#     5 = CIGAR	String	CIGAR String
#     6 = RNEXT	String	Ref. name of the mate/next read
#     7 = PNEXT	Int	Position of the mate/next read
#     8 = TLEN	Int	observed Template LENgth
#     9 = SEQ	String	segment SEQuence
#     10 = QUAL	String	ASCII of Phred-scaled base QUALity+33
#
# If -u is present, L1R1 and L2R1 are checked. If they are not in the UMI dictionary, that line is written to unknown.sam
#
# Use match_chromosome(L1R1[2], L1R2[2]), if this function returns False, write L1R1 to samfile_deduped.sam. The line after the original L1R1 becomes L1R1, and the algorithm restarts from Deduping. If it returns True, continue.
#
# Use match_UMI(L1R1[0], L1R2[0]), if they return False, iterate L2R1 to next line down, and restart from Deduping. If it returns True, continue.
#
# Use match_strand(L1R1[1], L2R1[1]), if this function returns False, iterate L2R1 to next line down, and restart from Deduping. If it returns True, continue.
#
# Use match_position(L1R1[2], L1R2[2]), if this function returns False, iterate L2R1 to next line down, and restart from Deduping. If it returns True, continue.
#
# If all functions have returned True, compare QSCORES (L1R1[4] and L2R1[4]), and keep the read with the highest number. If they are equivalent, discard L2R1. Write discarded entry to duplicates.txt if that option was selected. Continue from Deduping with line that was kept as L1R1.
#
# After reaching the end of the SAM file, close SAM file and any open output files.
