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
				_,rid,_,_,loc,seq,*_ = eval(line)
				values = products.get(rid, None)
				if values:
					locus = Locus(rid, seq)
					# fix negative numbers in the genbank feature locations
					#if '-1' in join:
					#	join = join.replace('-1', str(locus.length()))
					#feature = locus.read_feature('CDS ' + join)
					a,d = loc.split('..')
					b,c = values[6].split(':')[0].split('..')
					pairs = [[a,b],[c,d]]
					feature = locus.add_feature('CDS', +1, pairs)
					feature.tags['recode'] = ''
					outfile = open(rid, 'w')
					locus.write(outfile)
					outfile.close()
	
			elif line.startswith('INSERT INTO products VALUES '):
				line = line.rstrip().replace('INSERT INTO products VALUES ','')[:-1]
				values = eval(line)
				if values[6] and 'ribosomal_frameshift' in values[6]:
					products[values[3]] = values
	
