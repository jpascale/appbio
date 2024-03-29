#!/usr/bin/env python2
from Bio import SeqIO
import sys
import math
import time
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

WINDOW_DISCARD = 0
WINDOW_SEQ_TRIMMED = 1
WINDOW_RETURN_WHOLE_SEQ = 2

def boxplot(data, x=0):
	#import pdb; pdb.set_trace()
	sorted_data = np.array(data.items())
	#sorted_data = np.sort(sorted_data, 0)
	values = sorted_data[:,0]
	freqs = sorted_data[:,1]
	freqs = np.cumsum(freqs)
	freqs = freqs*1./np.max(freqs)

	#get 25%, 50%, 75% percentiles
	idx = np.searchsorted(freqs, [0.25, 0.5, 0.75])
	p25, p50, p75 = values[idx]
	vmin, vmax = values.min(), values.max()

	ax = plt.gca()
	l,r = -0.2+x, 0.2+x
	#plot boxes
	plt.plot([l,r], [p50, p50], 'k')
	plt.plot([l, r, r, l, l], [p25, p25, p75, p75, p25], 'k')
	plt.plot([x,x], [p75, vmax], 'k')
	plt.plot([x,x], [p25, vmin], 'k')

	return vmin, p25, p50, p75, vmax

class StatsHolder(object):

	def __init__(self):
		self.discard_counter = 0
		self.accept_counter = 0
		self.bp = dict()

	def start_timer(self):
		self.start_time = time.time()

	def stop_timer(self):
		self.elapsed_time = time.time() - self.start_time

	def discard_sequence(self):
		self.discard_counter += 1

	def accept_sequence(self):
		self.accept_counter += 1

	def add_bp_metric(self, start, stop, lenn):
		percentage = int(round((lenn - (stop-start))/float(lenn) * 100))
		if percentage in self.bp:
			self.bp[percentage] += 1
		else:
			self.bp[percentage] = 1

	def print_stats(self):
		print "Elapsed time: " + str(self.elapsed_time) + " seconds."

		total_seq = self.accept_counter + self.discard_counter
		if total_seq > 0:
			print "Discarded sequences: " + str(self.discard_counter) + "/" + str(total_seq) + " (" + str(self.discard_counter/float(total_seq) * 100) + "%)"
		sys.stdout.write('\n')
		print "Discarded bases distribution"

		vmin, p25, p50, p75, vmax = boxplot(self.bp,0)
		#import pdb; pdb.set_trace()
		print "vmin: " + str(vmin)
		print "q25: " + str(p25)
		print "q50: " + str(p50)
		print "q75: " + str(p75)
		print "vmax: " + str(vmax)
		plt.xlim(-0.5,1.5)
		plt.show()

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

def window_algorithm(record):
	"""
	Sildes a window of length = window_length across sequence and calculate
	the mean value within that window. If the mean value drops below THRESHOLD
	at a certain point, then a record will be returned representing the sequence
	with the 3' end removed at that point. If the trimmed read length is lower
	than READ_LENGTH_DEFAULT, the read will be discarded.
	"""
	seq = record.seq

	#Check if read is shorter than read length
	if len(seq) < READ_LENGTH_DEFAULT:
		return None, None, True
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
				pos = i + window_len/2
			break

		elif mean >= THRESHOLD:
			i += 1
		else:
			if i < READ_LENGTH_DEFAULT:
				break
			else:
				sub_rec = WINDOW_SEQ_TRIMMED
				pos = i + window_len/2
				break
   
	return sub_rec, pos, False

#def mott_algorithm_old(record):
#    """
#    Iterates over the sequence incorporating nucleotides with quality above the threshold.
#    Those nucleotides under the threshold surrounded by high quality nucleotides are 
#    incorporated with a tolerance of TOLERANCE_LEN. If the TOLERANCE_LEN is reached,
#    the subsequence is a candidate to be returned if there is no other candidate that is 
#    longer. If the trimmed length is  lower than READ_LENGTH_DEFAULT, the read will be discarded
#    """


def cumsum_algorithm(record):
	seq = record.seq
	#Check if read is horter than read length 
	if len(seq) < READ_LENGTH_DEFAULT:
		return None, None, True
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
		return 0, 0, False

	highest_p_pos = start_pos = i

	#Find the peak
	while i < len(seq):
		delta_p += threshold_prob - get_prob_by_quality(qual[i])
		if delta_p >= highest_p:
			highest_p = delta_p
			highest_p_pos = i
		i += 1
	
	if highest_p_pos - start_pos + 1 > READ_LENGTH_DEFAULT:
		return start_pos, highest_p_pos, None
	else:
		return 0, 0, False


def get_window_mean(qual):
	"""
	Calculates window mean.
	"""
	return float(sum(qual))/len(qual)

def get_prob_by_quality(quality):
	''' q = -10 * log10(p) '''
	return math.pow(10, quality / -10.0)

def main():
	"""
	Parses through fastq file and trims each record using the sliding window algorithm.
	"""
	setting_variables() #Obtain global variables from user input

	if os.stat(filename).st_size == 0:
		sys.exit("file is empty")
	if STATS:
		stats = StatsHolder()
		stats.start_timer()

	short_reads = 0 #Counter for reads shorter than default read length

	try:
		with open(filename) as ih, open(outputfile, 'w') as oh:
			if ALGORITHM == "cumsum":
				sys.stderr.write("Running CumSum algorithm\n")
				for record in SeqIO.parse(ih, 'fastq'):
					base, i, short = cumsum_algorithm(record)
					if short:
							short_reads += 1
							continue
					if not i == 0:
						oh.write(record[base:i+1].format('fastq'))

						if STATS:
							stats.accept_sequence()
							stats.add_bp_metric(base, i, len(record))

					else:
						if STATS:
							stats.discard_sequence()
			else:
				sys.stderr.write("Running sliding window algorithm\n")
				for record in SeqIO.parse(ih, 'fastq'):
					trimmed_seq, pos, short = window_algorithm(record)
					if short:
							short_reads += 1
							continue
					if trimmed_seq == WINDOW_SEQ_TRIMMED:
						oh.write(record[0:pos].format('fastq'))
						
						if STATS:
							stats.accept_sequence()
							stats.add_bp_metric(0, pos, len(record))


					elif trimmed_seq == WINDOW_RETURN_WHOLE_SEQ:
						oh.write(record.format('fastq'))
						
						if STATS:
							stats.accept_sequence()
							stats.add_bp_metric(0, len(record), len(record))

					else:
						if STATS:
							stats.discard_sequence()

	except NameError:      # If no output filename is specified, write to stdout.
		with open(filename) as ih: 
			try:
				if ALGORITHM == "cumsum":
					sys.stderr.write("Running CumSum algorithm\n")

					for record in SeqIO.parse(ih, 'fastq'):
						base, i, short = cumsum_algorithm(record)
						if short:
							short_reads += 1
							continue
						if not i == 0:
							sys.stdout.write(record[base:i+1].format('fastq'))

							if STATS:
								stats.accept_sequence()
								stats.add_bp_metric(base, i, len(record))

						else:
							if STATS:
								stats.discard_sequence()

				else:
					sys.stderr.write("Running sliding window algorithm\n")
					for record in SeqIO.parse(ih, 'fastq'):
						trimmed_seq, pos, short = window_algorithm(record)
						if short:
							short_reads += 1
							continue
						elif trimmed_seq == WINDOW_SEQ_TRIMMED:
							sys.stdout.write(record[0:pos].format('fastq'))

							if STATS:
								stats.accept_sequence()
								stats.add_bp_metric(0, pos, len(record))

						elif trimmed_seq == WINDOW_RETURN_WHOLE_SEQ:
							sys.stdout.write(record.format('fastq'))

							if STATS:
								stats.accept_sequence()
								stats.add_bp_metric(0, len(record), len(record))	
						else:
							if STATS:		
								stats.discard_sequence()
			except ValueError as e: #Obtain raised error from Biopython
				print(e)
				print("Please check that the file is correctly formatted.")
				sys.exit()

	if STATS:
		stats.stop_timer()
		stats.print_stats()
	if short_reads != 0:
		sys.stderr.write("Warning: " + str(short_reads) + " read(s) shorter than default read length\n")

if __name__ == "__main__":
	main()