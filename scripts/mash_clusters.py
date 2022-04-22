import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from itertools import chain as chain
from math import log
import collections
import math

def is_valid_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


usage = '%s [-opt1, [-opt2, ...]] file1 file2' % __file__
parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
parser.add_argument('infile', type=is_valid_file, help='input file')
parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
parser.add_argument('-c', '--cutoff', help='The minimum cutoff similarity [0.95]', type=float, default=0.95)
args = parser.parse_args()


names = dict()
clusters = list()
where = dict()

with open(args.infile) as fp:
	next(fp)
	for read2, line in enumerate(fp, start=0):
		name, *sims = line.rstrip().split('\t')
		names[read2] = name
		for read1,s in enumerate(sims, start=0):
			if float(s) < 1 - args.cutoff:
				#print(names[read1], names[read2], s)
				if read1 in where and read2 in where:
					if where[read1] != where[read2]:
						clusters.remove(where[read2])
						for item in where[read2]:
							where[read1].append(item)
							where[item] = where[read1]
				elif read1 in where:
					where[read1].append(read2)
					where[read2] = where[read1]
				elif read2 in where:
					where[read2].append(read1)
					where[read1] = where[read2]
				else:
					lis = [read1, read2]
					clusters.append(lis)
					where[read1] = lis
					where[read2] = lis

for cluster in clusters:
	print(len(cluster), end='\t')
	for item in sorted(cluster, key = lambda item : names[item] ):
		print(names[item], end='\t')
		del names[item]
	print()

for item in names:
	print(1, names[item], sep='\t')







	
