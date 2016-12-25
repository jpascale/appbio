#!/usr/bin/env python2
from Bio import SeqIO
import sys
#import time
import argparse

WINDOW_SEQ_TRIMMED = 1
WINDOW_RETURN_WHOLE_SEQ = 2

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
    parser.add_argument("-tl", "--tolerance", type=int,
                        help="Set tolerance error len for the mott algorithm. Default is 4.",
                        action="store")
    parser.add_argument("-a", "--algorithm", type=str,
                        help="Use \'mott\' for Mott algorithm usage. Otherwise window algorithm will be used.",
                        action="store")
    args = parser.parse_args()
    
    global filename
    filename = args.input
    
    global THRESHOLD   
    if args.threshold:
        THRESHOLD = args.threshold
    else:
        THRESHOLD = 20

    global TOLERANCE_LEN
    if args.tolerance:
        TOLERANCE_LEN = args.tolerance
    else:
        TOLERANCE_LEN = 6
        
    global ALGORITHM
    ALGORITHM = "mott" if args.algorithm == "mott" else "window"

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
                sub_rec = WINDOW_RETURN_WHOLE_SEQ
                break
        elif mean >= THRESHOLD:
            i += 1
        else:
            if i < READ_LENGTH_DEFAULT:
                break
            else:
                sub_rec = WINDOW_SEQ_TRIMMED
                pos = i
                break
   
    return sub_rec, pos

def mott_algorithm(record):
    """
    Iterates over the sequence incorporating nucleotides with quality above the threshold.
    Those nucleotides under the threshold surrounded by high quality nucleotides are 
    incorporated with a tolerance of TOLERANCE_LEN. If the TOLERANCE_LEN is reached,
    the subsequence is a candidate to be returned if there is no other candidate that is 
    longer.
    """

    seq = record.seq
    qual = record.letter_annotations["phred_quality"]

    err_count = 0

    base = i = 0
    max_seq_base = max_seq_i = 0 #Keeps record of the longest subsequence

    while i < len(seq):

        if qual[i] >= THRESHOLD:
            err_count = 0 #Reset quality error counting.
        else:
            err_count += 1
            if err_count == TOLERANCE_LEN and base + i - TOLERANCE_LEN > max_seq_i - max_seq_base: ## Check if we have a new candidate
                    max_seq_base = base
                    max_seq_i = i - TOLERANCE_LEN + 1
                    base = i

        i += 1
    print "DEBUG: _________NEW_________________"
    print seq
    print "*****"

    if max_seq_base == 0 and max_seq_i == 0: ## The whole sequence is returned
        print seq
        return 0, len(seq)
    elif i - base > max_seq_i - max_seq_base:
        print seq[base:i]
        return base, i 
    else:
        print seq[max_seq_base:max_seq_i]
        return max_seq_base, max_seq_i



def get_window_mean(qual):
    """
    Calculates window mean.
    """
    return float(sum(qual))/len(qual)

def main():
    #start = time.time()
    """
    Parses through fastq file and trims each record using the sliding window algorithm.
    """
    setting_variables()
    try:
        with open(filename) as ih, open(outputfile, 'w') as oh:
            if ALGORITHM == "mott":
                print "Running mott"
                for record in SeqIO.parse(ih, 'fastq'):
                    base, i = mott_algorithm(record)
                    oh.write(record[base:i].format('fastq'))
            else:
                for record in SeqIO.parse(ih, 'fastq'):
                    trimmed_seq, pos = window_algorithm(record)
                    if trimmed_seq == WINDOW_SEQ_TRIMMED:
                        oh.write(record[0:pos].format('fastq'))
                    elif trimmed_seq == WINDOW_RETURN_WHOLE_SEQ:
                        oh.write(record.format('fastq'))
    except NameError:
        with open(filename) as ih: 
            if ALGORITHM == "mott":
                for record in SeqIO.parse(ih, 'fastq'):
                    base, i = mott_algorithm(record)
                    oh.write(record[base:i].format('fastq'))
            else:
                for record in SeqIO.parse(ih, 'fastq'):
                    trimmed_seq, pos = window_algorithm(record)
                    if trimmed_seq == WINDOW_SEQ_TRIMMED:
                        sys.stdout.write(record[0:pos].format('fastq'))
                    elif trimmed_seq == WINDOW_RETURN_WHOLE_SEQ:
                        sys.stdout.write(record.format('fastq'))
            
    #end = time.time()
    #print(end - start)
if __name__ == "__main__":
	main()