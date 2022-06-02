#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from types import MethodType
from termcolor import colored

from genbank.locus import Locus

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x
def _print(self, item):
	if isinstance(item, str):
		self.write(item)
	else:
		self.write(str(item))


if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-d', '--dump', action="store_true")
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	NULL = None
	products = dict()
	with open(sys.argv[1]) as fp:
		for line in fp:
			if line.startswith('INSERT INTO molecules VALUES'):
				line = line.rstrip().replace('INSERT INTO molecules VALUES ','')[:-1]
				values = eval(line)
				rid = values[1]
				if rid in products:
					fs = products[rid]
					if fs and 'ribosomal_frameshift' in fs:
						locus = Locus(rid, values[5])
						#locus.read_feature('CDS ' + join)
						a0,a1 = fs.split(':')[0].split('..')
						pairs = [['1',str(int(a0))],[str(int(a1)-3),str(locus.length())]]
						feature = locus.add_feature('CDS', +1, pairs)
						feature.tags['recode'] = fs
						outfile = open(rid, 'w')
						locus.write(outfile)
						outfile.close()
	
			elif line.startswith('INSERT INTO products VALUES '):
				line = line.rstrip().replace('INSERT INTO products VALUES ','')[:-1]
				values = eval(line)
				products[values[3]] = values[6]
	
