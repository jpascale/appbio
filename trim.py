#!/usr/bin/env python2
from Bio import SeqIO
import sys
#import time
import argparse

def setting_variables():
    """
    Define parameters from user settings.
    """
    parser = argparse.ArgumentParser()
    #Positional argument
    parser.add_argument("input", 
                        help="Give path to input file. Only fastq format allowed.")    
    #Optional arguments
    parser.add_argument("-t", "--threshold", type=int,
                        help="Set quality threshold for trimming.",
                        action="store")
    parser.add_argument("-l", "--length", type=int,
                        help="Set minimum read length. Any reads smaller than this value will be discarded.",
                        action="store")
    parser.add_argument("-w", "--window", type=int,
                        help="Set window width, i.e. the number of bases to include in each window.",
                        action="store")
    parser.add_argument("-p", "--percentage_window", type=float,
                        help="""Set window size as percentage of total read instead of fixed number of base pairs. 
                        Minimum window width is set to 3 base pairs. If a default window size is given with -w, 
                        the largest window size will be selected.""")
    parser.add_argument("-out", "--outputname", type=str,
                        help="output file name",
                        action="store")
    args = parser.parse_args()
    
    global filename
    filename = args.input
    
    global THRESHOLD   
    if args.threshold:
        THRESHOLD = args.threshold
    else:
        THRESHOLD = 20
        
    global READ_LENGTH_DEFAULT
    if args.length:
        READ_LENGTH_DEFAULT = args.length
    else:
        READ_LENGTH_DEFAULT = 20
        
    global WINDOW_DEFAULT
    if args.window:
        WINDOW_DEFAULT = args.window
    else:
        WINDOW_DEFAULT = 3
        
    global WINDOW_PERCENT
    if args.percentage_window:
        WINDOW_PERCENT = args.percentage_window
    else:
        WINDOW_PERCENT = 0.0 #Default percentage set to 0 so that window width provded by user will be constant
    
    if args.outputname:
        global outputfile
        outputfile = args.outputname

#TODO: Test
def window_algorithm(record):
    """
    Sildes a window of length = window_length across sequence and calculate
    the mean value within that window. If the mean value drops below THRESHOLD
    at a certain point, then a record will be returned representing the sequence
    with the 3' end removed at that point. If the trimmed read length is lower
    than READ_LENGTH_DEFAULT, the read will be discarded.
    """
    seq = record.seq
    qual = record.letter_annotations["phred_quality"]
    
    window_len = max(WINDOW_DEFAULT, int(round(len(seq) * WINDOW_PERCENT)))
    sub_rec = 0
    pos = 0
    i = 0
    
    while i + window_len <= len(seq):
        subseq_qual = qual[i:i+window_len]      #Obtain quality values in window
        mean = get_window_mean(subseq_qual)     #Calculate mean value in window
        if i + window_len == len(seq):
                sub_rec = 2
                break
        elif mean >= THRESHOLD:
            i += 1
        else:
            if i < READ_LENGTH_DEFAULT:
                break
            else:
                sub_rec = 1
                pos = i
                break
   
    return sub_rec, pos

def get_window_mean(qual):
    """
    Calculates window mean.
    """
    return float(sum(qual))/len(qual)

def main():
    #start = time.time()
    """
    Parses trhough fastq file and trims each record using the sliding window algorithm.
    """
    setting_variables()
    try:
        with open(filename) as ih, open(outputfile, 'w') as oh:
            for record in SeqIO.parse(ih, 'fastq'):
                trimmed_seq, pos = window_algorithm(record)
                if trimmed_seq == 1:
                    oh.write(record[0:pos].format('fastq'))
                elif trimmed_seq == 2:
                    oh.write(record.format('fastq'))
    except NameError:
        with open(filename) as ih:
            for record in SeqIO.parse(ih, 'fastq'):
                trimmed_seq, pos = window_algorithm(record)
                if trimmed_seq == 1:
                    sys.stdout.write(record[0:pos].format('fastq'))
                elif trimmed_seq == 2:
                    sys.stdout.write(record.format('fastq'))
            
    #end = time.time()
    #print(end - start)
if __name__ == "__main__":
	main()