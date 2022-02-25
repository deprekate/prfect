import sys
#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 


def get_stop(seq):
	for i in range(0,len(seq)-2,3):
		codon = seq[i:i+3]
		if codon in ['TAA','TAG','TGA']:
			return i
	return None
		


with open(sys.argv[1]) as fp:
	for line in fp:
		cols = line.rstrip().split('\t')
		if int(cols[4]) / int(cols[3]) < 0.9:
			continue
		seq = cols[13].upper()
		print(seq, get_stop(seq))

