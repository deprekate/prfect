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
from packaging import version

sys.path.pop(0)

from genbank.feature import Feature
from prfect.file import File
import pandas as pd

# sklearn and model persisitence is iffy
import sklearn
if version.parse(sklearn.__version__) < version.parse('1.0.0'):
	from sklearn.experimental import enable_hist_gradient_boosting
	path = pkg_resources.resource_filename('prfect', 'clf.0.24.0.pkl')
elif version.parse(sklearn.__version__) < version.parse('1.1.0'):
	path = pkg_resources.resource_filename('prfect', 'clf.1.0.pkl')
elif version.parse(sklearn.__version__) < version.parse('1.1.1'):
	path = pkg_resources.resource_filename('prfect', 'clf.1.1.0.pkl')
else:
	path = pkg_resources.resource_filename('prfect', 'clf.1.1.1.pkl')
from sklearn.ensemble import HistGradientBoostingClassifier
clf = pickle.load(open(path, 'rb'))

def strr(x):
    if isinstance(x, float):
        return str(round(x,4))
    else:
        return str(x)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x


def alert(args, label, last, curr, metrics):
	sys.stderr.write(colored("ribo frameshift detected in " + args.infile + "\n", 'red') )
	args.outfile.print("\n")
	args.outfile.print("     CDS             join(%s..%s,%s..%s)" % (last.left(), last.right(), curr.left(), curr.right()))
	args.outfile.print("\n")
	args.outfile.print("                     /ribosomal_slippage=%s" % metrics['DIR']  )
	args.outfile.print("\n")
	args.outfile.print("                     /slippery_sequence=%s" % metrics['W'] + metrics['E0'] + metrics['P0'] + metrics['A0'] )
	args.outfile.print("\n")
	args.outfile.print("                     /motif=%s" % args.locus.number_motif(metrics['MOTIF']).__name__  )
	args.outfile.print("\n")
	args.outfile.print("                     /label=%s" % label )
	args.outfile.print("\n")
	if 'product' in last.tags or 'product' in curr.tags:
		args.outfile.print("                     /product=%s,%s" % (last.tags.get('product',''),curr.tags.get('product','')) )
		args.outfile.print("\n")

flag = True
def dump(args, label, last, curr, metrics):
	# this is to set only frameshifts that occur within 10bp
	if 10 > ((last.right() + curr.left()) / 2 - metrics['LOC']):
		metrics['LABEL'] = 1

	global flag
	if flag:
		args.outfile.print('GENOME\t')
		args.outfile.print('\t'.join(map(str,metrics.keys())))
		args.outfile.print('\n')
		flag = False
	args.outfile.print(args.locus.name)
	args.outfile.print('\t')
	args.outfile.print('\t'.join(map(strr,metrics.values())))
	args.outfile.print('\n')

def _print(self, item):
	if isinstance(item, str):
		self.write(item)
	else:
		self.write(str(item))

def has_prf(metrics):
	global clf
	row = pd.DataFrame.from_dict(metrics,orient='index').T
	if clf.predict(row.loc[:,clf.feature_names_in_])[0] == metrics['DIR']:
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
		args.locus = locus
		_last = _curr = None
		for feature in locus:
			#if feature.is_type('CDS') and feature.is_joined() and '1' not in sum(feature.pairs, ()) and len(feature.pairs)==2 and int(feature.pairs[1][0])-int(feature.pairs[0][1]) < 100:
			if feature.is_type('CDS') and feature.is_joined() and len(feature.pairs)==2 and abs(int(feature.pairs[1][0])-int(feature.pairs[0][1])) < 10:
				#sys.stderr.write(colored("Genome already has a joined feature:\n", 'red') )
				#feature.write(sys.stderr)
				#sys.stderr.write(colored("...splitting the feature into two for testing\n\n", 'red') )
				_last = Feature(feature.type, feature.strand, [feature.pairs[0]], locus, feature.tags)
				_curr = Feature(feature.type, feature.strand, [feature.pairs[1]], locus, feature.tags)
				for metrics in locus.get_metrics(_last, _curr):
					if args.dump:
						dump(args, 1, _last, _curr, metrics)
					elif has_prf(metrics):
						alert(args, 1, _last, _curr, metrics)
				_last = None
			elif feature.is_type('CDS') and len(feature.pairs)==1:
				#continue
				if _last and _last.strand==feature.strand:
					for metrics in locus.get_metrics(_last, feature):
						if args.dump:
							dump(args, 0, _last, feature, metrics)
						elif has_prf(metrics):
							alert(args, 0, _last, feature, metrics)
				_last = feature
	

