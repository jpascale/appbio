#!/usr/bin/env python2
from Bio import SeqIO
import sys
import time

THRESHOLD = 20              #Quality threshold
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
    print(record.__hash__)
    
    window_len = max(WINDOW_DEFAULT, int(round(len(seq) * WINDOW_PERCENT)))
    sub_rec = 0
    rem_rec = 0
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
                sub_rec = 1
                break
            else:
                sub_rec = record[0:i]
                rem_rec = record[i:]
                break

   
    return sub_rec, rem_rec#, discarded_reads, removed_reads

def get_window_mean(qual):
    """
    Calculates window mean.
    """
    return float(sum(qual))/len(qual)

def main():
    start = time.time()
    """
    Parses trhough fastq file and trims each record using the sliding window algorithm.
    """
    filename = sys.argv[1]
    output = sys.argv[2]
    with open(filename) as ih, open(output, 'w') as oh, open('discarded', 'w') as df, open('removed', 'w') as rf:
        for record in SeqIO.parse(ih, 'fastq'):
            trimmed_seq, removed = window_algorithm(record)

            if type(trimmed_seq) == type(record):
                oh.write(trimmed_seq.format('fastq'))
                rf.write(removed.format('fastq'))
            elif trimmed_seq == 1:
                df.write(record.format('fastq'))
            elif trimmed_seq == 2:
                oh.write(record.format('fastq'))
            
    end = time.time()
    print(end - start)
if __name__ == "__main__":
	main()