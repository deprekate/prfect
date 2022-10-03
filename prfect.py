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
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6
import pandas as pd
import numpy as np

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
clf = None

def strr(x):
    if isinstance(x, float):
        return str(round(x,5))
    else:
        return str(x)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x


def alert(args, last, curr, metrics):
	# this is to set only frameshifts that occur within 10bp
	#if label and 10 > ((last.right() + curr.left()) / 2 - metrics['LOC']):

	pairs = [[last.left(), last.right()], [curr.left(), curr.right()]]
	if last.strand > 0:
		pairs[0][-1] = metrics['LOC'] + 2
		pairs[-1][0] = metrics['LOC'] + 3 + metrics['DIR']
	pairs = [list(map(str, lis)) for lis in pairs] 
	feature = Feature('CDS', curr.strand, pairs, args.locus)
	
	feature.tags['ribosomal_slippage'] =  metrics['DIR']
	feature.tags['motif'] = args.locus.number_motif(metrics['MOTIF']).__name__
	feature.tags['label'] = metrics['LABEL']
	feature.tags['locus'] = args.locus.name()
	feature.tags['location'] = metrics['LOC']
	feature.tags['translation'] = feature.translation()
	if 'product' in last.tags or 'product' in curr.tags:
		feature.tags['product'] = '"' + last.tags.get('product','').replace('"','') + ';' + curr.tags.get('product','').replace('"', '') + '"'
	if args.format == 'feature':
		feature.write(args.outfile)

flag = True
def dump(args, last, curr, metrics):
	global flag
	if flag:
		args.outfile.print('GENOME\t')
		args.outfile.print('\t'.join(map(str,metrics.keys())))
		args.outfile.print('\n')
		flag = False
	args.outfile.print(args.locus.name())
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
	prob = clf.predict_proba(row.loc[:,clf.feature_names_in_])
	metrics['pred'] = clf.classes_[np.argmax(prob)]
	metrics['prob'] = np.max(prob)
	if metrics['pred'] == metrics['DIR']:
		return True
		

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-d', '--dump', action="store_true")
	parser.add_argument('-p', '--param', type=str) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	parser.add_argument('-m', '--model', type=str) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='feature', choices=['tabular','genbank','feature'])
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)
	
	if args.dump:
		pass
	elif args.model and not args.dump:
		clf = pickle.load(open(args.model, 'rb'))
	else:
		clf = pickle.load(open(path, 'rb'))

	genbank = File(args.infile)
	for name,locus in genbank.items():
		#for codon,rarity in locus.codon_rarity().items():print(codon, rarity, sep='\t')
		locus.init(args)
		args.locus = locus
		locus.args = args
		_last = _curr = None
		for feature in locus:
			best = dict()
			if feature.is_type('CDS') and feature.is_joined() and len(feature.pairs)==2 and abs(int(feature.pairs[1][0])-int(feature.pairs[0][1])) < 10:
				#sys.stderr.write(colored("Genome already has a joined feature:\n", 'red') )
				#feature.write(sys.stderr)
				#sys.stderr.write(colored("...splitting the feature into two for testing\n\n", 'red') )
				_last = Feature(feature.type, feature.strand, [feature.pairs[0]], locus, feature.tags)
				_curr = Feature(feature.type, feature.strand, [feature.pairs[1]], locus, feature.tags)
				for metrics in locus.get_metrics(_last, _curr):
					#metrics['LABEL'] = 1  if 10 > abs((_last.right() + _curr.left()) / 2 - metrics['LOC']) else 0
					metrics['LABEL'] = 1
					if 10 < abs((_last.right() + _curr.left()) / 2 - metrics['LOC']):
						metrics['LABEL'] = 0
						#continue
					if args.dump:
						dump(args, _last, _curr, metrics)
					elif has_prf(metrics):
						if not best or metrics['prob'] > best['prob']:
							best = metrics
				if best:
					alert(args, _last, _curr, best)
				_last = None
			elif feature.is_type('CDS') and feature.is_joined() and len(feature.pairs)==3:
				for pair1, pair2 in zip(feature.pairs, feature.pairs[1:]):
					_last = Feature(feature.type, feature.strand, [pair1], locus, feature.tags)
					_curr = Feature(feature.type, feature.strand, [pair2], locus, feature.tags)
					for metrics in locus.get_metrics(_last, _curr):
						metrics['LABEL'] = 1
						if args.dump:
							dump(args, _last, _curr, metrics)
						elif has_prf(metrics):
							if not best or metrics['prob'] > best['prob']:
								best = metrics
					if best:
						alert(args, _last, _curr, best)
				_last = None
			elif feature.is_type('CDS') and len(feature.pairs)==1:
				#continue
				if _last and _last.strand==feature.strand:
					for metrics in locus.get_metrics(_last, feature):
						if args.dump:
							dump(args, _last, feature, metrics)
						elif has_prf(metrics):
							if not best or metrics['prob'] > best['prob']:
								best = metrics
					if best:
						alert(args, _last, feature, best)
				_last = feature
			if not best and not args.dump:
				feature.write(args.outfile)
				pass
	

