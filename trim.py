#!/usr/bin/env python2
from Bio import SeqIO
import sys
import math
import time
import argparse

WINDOW_DISCARD = 0
WINDOW_SEQ_TRIMMED = 1
WINDOW_RETURN_WHOLE_SEQ = 2

class StatsHolder(object):

    def __init__(self):
        self.discard_counter = 0
        self.accept_counter = 0

    def start_timer(self):
        self.start_time = time.time()

    def stop_timer(self):
        self.elapsed_time = time.time() - self.start_time

    def discard_sequence(self):
        self.discard_counter += 1

    def accept_sequence(self):
        self.accept_counter += 1

    def print_stats(self):
        print "Elapsed time: " + str(self.elapsed_time) + " seconds."

        total_seq = self.accept_counter + self.discard_counter
        if total_seq > 0:
            print "Discarded sequences: " + str(self.discard_counter) + "/" + str(total_seq) + " (" + str(self.discard_counter/float(total_seq) * 100) + "%)"



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
    parser.add_argument("-a", "--algorithm", type=str,
                        help="""Use \'cumsum\' for Cumsum algorithm usage.
                        Parameters with cumsum algorithm: Threshold (-t int), Minimum read lenght (-l int), Output filename (-out str). 
                        Otherwise sliding window algorithm will be used.""",
                        action="store")
    parser.add_argument("-s", "--stats",
                        help="See statistics",
                        action="store_true")
    args = parser.parse_args()
    
    global filename
    filename = args.input

    global ALGORITHM
    ALGORITHM = "cumsum" if args.algorithm == "cumsum" else "window"
    
    global THRESHOLD   
    if ALGORITHM == "cumsum":
        THRESHOLD = 20 # Threshold for the cumsum
    else:
        THRESHOLD = 20 # Threshold for the window slide

    if args.threshold:
        THRESHOLD = args.threshold

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

    global STATS
    if args.stats:
        STATS = True
    else:
        STATS = False

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
    sub_rec = WINDOW_DISCARD
    pos = 0
    i = 0
    
    while i + window_len <= len(seq):
        subseq_qual = qual[i:i + window_len]      #Obtain quality values in window
        mean = get_window_mean(subseq_qual)     #Calculate mean value in window
        if i + window_len == len(seq):

            # Patch that solves the problem with the final part of the sequence
            if mean >= THRESHOLD:
                #Patch: Only add entire sequence if last nucleotides are not under the quality
                if qual[-1] < THRESHOLD:
                    j = -1
                    while qual[j] < THRESHOLD:
                        j -= 1
                    pos = len(seq) + j
                    sub_rec = WINDOW_SEQ_TRIMMED
                else:
                    sub_rec = WINDOW_RETURN_WHOLE_SEQ

            else:
                sub_rec = WINDOW_SEQ_TRIMMED
                pos = i
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

#def mott_algorithm_old(record):
#    """
#    Iterates over the sequence incorporating nucleotides with quality above the threshold.
#    Those nucleotides under the threshold surrounded by high quality nucleotides are 
#    incorporated with a tolerance of TOLERANCE_LEN. If the TOLERANCE_LEN is reached,
#    the subsequence is a candidate to be returned if there is no other candidate that is 
#    longer. If the trimmed length is  lower than READ_LENGTH_DEFAULT, the read will be discarded
#    """

#TODO: document
def cumsum_algorithm(record):
    seq = record.seq
    qual = record.letter_annotations["phred_quality"]

    highest_p = 0
    highest_p_pos = 0
    delta_p = 0

    threshold_prob = get_prob_by_quality(THRESHOLD)

    # Remove sequences from the beginning
    i = 0
    while i < len(seq) and get_prob_by_quality(qual[i]) > threshold_prob:
        i += 1

    # Limit case in which all the nucleotides are under the threshold, seq is discarded
    if i == len(seq):
        return 0, 0

    highest_p_pos = start_pos = i

    #Find the peak
    while i < len(seq):
        delta_p += threshold_prob - get_prob_by_quality(qual[i])
        if delta_p > highest_p:
            highest_p = delta_p
            highest_p_pos = i
        i += 1
    
    if highest_p_pos - start_pos + 1 > READ_LENGTH_DEFAULT:
        print str(qual[start_pos:start_pos+10]) + " " + str(start_pos) + " " + str(highest_p_pos) + " " + str(len(seq))
        return start_pos, highest_p_pos
    else:
        return 0, 0


def get_window_mean(qual):
    """
    Calculates window mean.
    """
    return float(sum(qual))/len(qual)

def get_prob_by_quality(quality):
    ''' q = -10 * log10(p) '''
    return math.pow(10, quality / -10.0)

def main():
    #start = time.time()
    """
    Parses through fastq file and trims each record using the sliding window algorithm.
    """
    setting_variables()
    if STATS:
        stats = StatsHolder()
        stats.start_timer()
    try:
        with open(filename) as ih, open(outputfile, 'w') as oh:
            if ALGORITHM == "cumsum":
                print "Running CumSum algorithm"
                for record in SeqIO.parse(ih, 'fastq'):
                    base, i = cumsum_algorithm(record)
                    if not i == 0:
                        stats.accept_sequence()
                        oh.write(record[base:i].format('fastq'))
                    else:
                        stats.discard_sequence()
            else:
                print "Running sliding window algorithm"
                for record in SeqIO.parse(ih, 'fastq'):
                    trimmed_seq, pos = window_algorithm(record)
                    if trimmed_seq == WINDOW_SEQ_TRIMMED:
                        oh.write(record[0:pos].format('fastq'))
                        stats.accept_sequence()
                    elif trimmed_seq == WINDOW_RETURN_WHOLE_SEQ:
                        oh.write(record.format('fastq'))
                        stats.accept_sequence()
                    else:
                        stats.discard_sequence()

    except NameError:
        print "NAME ERROR  NAME ERROR  NAME ERROR  NAME ERROR  NAME ERROR  NAME ERROR  NAME ERROR  NAME ERROR"
        with open(filename) as ih: 
            if ALGORITHM == "cumsum":
                print "Running CumSum algorithm"
                for record in SeqIO.parse(ih, 'fastq'):
                    base, i = cumsum_algorithm(record)
                    if not i == 0:
                        stats.accept_sequence()
                        sys.stdout.write(record[base:i].format('fastq'))
                    else:
                        stats.discard_sequence()

            else:
                print "Running sliding window algorithm"
                for record in SeqIO.parse(ih, 'fastq'):
                    trimmed_seq, pos = window_algorithm(record)
                    if trimmed_seq == WINDOW_SEQ_TRIMMED:
                        sys.stdout.write(record[0:pos].format('fastq'))
                        stats.accept_sequence()
                    elif trimmed_seq == WINDOW_RETURN_WHOLE_SEQ:
                        sys.stdout.write(record.format('fastq'))
                        stats.accept_sequence()
                    else:
                        stats.discard_sequence()

    if STATS:
        stats.stop_timer()
        stats.print_stats()
            
    #end = time.time()
    #print(end - start)
if __name__ == "__main__":
	main()