import math
import pdb

THRESHOLD = 20
WINDOW_PERCENT = 0.05
WINDOW_DEFAULT = 3 #Value used for window length as a minimum

def get_error_prob(q):
	return math.pow(10, ord(q) / -10.0)

ERROR_PROB_TOLERANCE = get_error_prob(THRESHOLD)

class Entry(object):

	#Respectively line 1, 2, 3 and 4
	def __init__(self, name, seq, plus_line, qual):
		self.name = name
		self.seq = seq
		self.plus_line = plus_line
		self.qual = qual


#TODO: Test
def window_algorithm(entry):
	seq = entry.seq
	qual = entry.qual

	window_len = max(WINDOW_DEFAULT, int(round(len(seq) * WINDOW_PERCENT)))
	saved_seq = []
	saved_qual = []

	i = 0

	while i + window_len < len(seq):
		subseq_qual = qual[i:i+window_len]
		mean = get_window_mean(subseq_qual)

		if mean < ERROR_PROB_TOLERANCE:
			saved_seq.append(seq[i])
			sabed_qual.append(qual[i])
			i += 1
		else:
			i += window_len #discards the whole window

	return Entry(entry.name, ''.join(saved_seq), entry.plus_line, ''.join(saved_qual))

def get_window_mean(qual):
	arr = [get_error_prob(c) for c in qual]
	return float(sum(arr))/len(arr)

def parse_file(filename, output, algorithm=None):
	#pdb.set_trace()
	with open(filename, 'r') as fd:
		while True:
			line = fd.readline().strip()

			if line == '': #EOF
				break

			name = line
			seq = fd.readline().strip()
			plus_line = fd.readline().strip()
			qual = fd.readline().strip()

			entry = Entry(name, seq, plus_line, qual)

			#Run algorithm over entry

			#Write
			write_entry(output, entry)
			

def write_entry(filename, entry):
	fd2 = open(filename, 'a')
	fd2.write(entry.name + '\n')
	fd2.write(entry.seq + '\n')
	fd2.write(entry.plus_line + '\n')
	fd2.write(entry.qual + '\n')
	fd2.close()

if __name__ == "__main__":
	parse_file("miseq1.fq", "output")