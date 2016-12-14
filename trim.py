import math
import pdb

def get_error_prob(q):
	return math.pow(10, 20 / -10.0)

class Entry(object):

	#Respectively line 1, 2, 3 and 4
	def __init__(self, name, seq, plus_line, qual):
		self.name = name
		self.seq = seq
		self.plus_line = plus_line
		self.qual = qual


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