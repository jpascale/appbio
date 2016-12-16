#!/usr/bin/env python2
from Bio import SeqIO
import sys

THRESHOLD = 15              #Quality threshold
WINDOW_PERCENT = 0.05       #Used to calculate window length from total length
WINDOW_DEFAULT = 3          #Value used for window length as a minimum
READ_LENGTH_DEFAULT = 20    #Minumum read length


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
    i = 0
    
    while i + window_len < len(seq):
        subseq_qual = qual[i:i+window_len]      #Obtain quality values in window
        mean = get_window_mean(subseq_qual)     #Calculate mean value in window
        if mean >= THRESHOLD:
		  i += 1
        else:
            if i < READ_LENGTH_DEFAULT:
                break
            else:
                sub_rec = record[0:i]
                break

   
    return sub_rec#, discarded_reads, removed_reads
 
def get_window_mean(qual):
    """
    Calculates window mean.
    """
    return float(sum(qual))/len(qual)

def main():
    """
    Parses trhough fastq file and trims each record using the sliding window algorithm.
    """
    filename = sys.argv[1]
    output = sys.argv[2]
    with open(filename) as ih, open(output, 'a') as oh:
        for record in SeqIO.parse(ih, 'fastq'):
            trimmed_seq = window_algorithm(record)

            if type(trimmed_seq) == type(record):
                oh.write(trimmed_seq.format('fastq'))
            else:
                continue
            

if __name__ == "__main__":
	main()