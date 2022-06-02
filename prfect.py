#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from subprocess import Popen, PIPE, STDOUT
from types import MethodType
from termcolor import colored
import pickle
import pkgutil
import pkg_resources

import pandas as pd

sys.path.pop(0)
from genbank.feature import Feature
#import prfect
from prfect.file import File
#from prfect.motif import Motif

#def extra(self, value=None):
#	return 'extra'
#locus.rare_codons = MethodType(rare_codons, locus)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def strr(x):
    if isinstance(x, float):
        return str(round(x,5))
    else:
        return str(x)

def alert(args, label, last, curr, features):
	sys.stderr.write(colored("ribo frameshift detected in " + args.infile + "\n", 'red') )
	args.outfile.print("\n")
	args.outfile.print("     CDS   join(%s..%s,%s..%s)" % (last.left(), last.right(), curr.left(), curr.right()))
	args.outfile.print("\n")
	args.outfile.print("           /product=%s" % last.tags['product'] )
	args.outfile.print("\n")
	args.outfile.print("           /product=%s" % curr.tags['product'] )
	args.outfile.print("\n")

flag = True
def dump(args, label, last, curr, features):
	global flag
	if flag:
		args.outfile.print('LABEL\tLASTL\tLASTR\tCURRL\tCURRR\t')
		args.outfile.print('\t'.join(map(str,features.keys())))
		flag = False
	args.outfile.print(label)
	args.outfile.print('\t')
	args.outfile.print(last.left())
	args.outfile.print('\t')
	args.outfile.print(last.right())
	args.outfile.print('\t')
	args.outfile.print(curr.left())
	args.outfile.print('\t')
	args.outfile.print(curr.right())
	args.outfile.print('\t')
	args.outfile.print('\t'.join(map(strr,features.values())))
	args.outfile.print('\n')

def _print(self, item):
	if isinstance(item, str):
		self.write(item)
	else:
		self.write(str(item))

path = pkg_resources.resource_filename('prfect', 'all.pkl')
#data = pkgutil.get_data(__name__, "all.pkl")
clf = pickle.load(open(path, 'rb'))
def has_prf(features):
	global clf
	row = pd.DataFrame.from_dict(features,orient='index').T
	row['RATIO'] = row['A1%'] / row['A0%']
	take = ['DIR', 'GC', 'N', 'RBS1','RBS2', 'A0%', 'A1%', 'RATIO', 'MOTIF', 'PROB', 'LF_35_6_RIGHT','HK_35_6_RIGHT','LF_40_6_RIGHT','HK_40_6_RIGHT']
	if clf.predict(row.loc[:,take])[0] == features['DIR']:
		return True
		

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-d', '--dump', action="store_true")
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	

	genbank = File(args.infile)
	for name,locus in genbank.items():
		#for codon,rarity in locus.codon_rarity().items():print(codon, rarity, sep='\t')
		locus.init()
		locus.args = args
		_last = _curr = None
		for feature in locus:
			if feature.is_type('CDS') and feature.is_joined() and '1' not in sum(feature.pairs, ()) and len(feature.pairs)==2 and int(feature.pairs[1][0])-int(feature.pairs[0][1]) < 100:
				a,b,c,d = map(int, sum(feature.pairs, () ))
				if abs(c-b) < 5:
					sys.stderr.write(colored("Genome already has a joined feature:\n", 'red') )
					feature.write(sys.stderr)
					sys.stderr.write(colored("...splitting the feature into two for testing\n\n", 'red') )
					b = b - (b-a-2)%3
					c = c + (d-c-2)%3
					a,b,c,d = map(str, [a,b,c,d])
					_last = Feature(feature.type, feature.strand, [[a,b]], locus, feature.tags)
					_curr = Feature(feature.type, feature.strand, [[c,d]], locus, feature.tags)
					for slip in locus.get_slips(_last, _curr):
						if args.dump:
							dump(args, 1, _last, _curr, slip)
						elif has_prf(slip):
							alert(args, 1, _last, _curr, slip)
					_last = None
			elif feature.is_type('CDS') and len(feature.pairs)==1:
				if _last and _last.strand==feature.strand:
					for slip in locus.get_slips(_last, feature):
						if args.dump:
							dump(args, 0, _last, feature, slip)
						elif False:
							alert(args, 0, _last, feature)
				_last = feature
	

