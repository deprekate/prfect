#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import os
import sys
import argparse
#from subprocess import Popen, PIPE, STDOUT
#from types import MethodType
#from termcolor import colored
from itertools import tee, chain
import pickle
import pkgutil
import pkg_resources
from packaging import version
import warnings
warnings.filterwarnings("ignore")
#warnings.filterwarnings("ignore", category=DeprecationWarning)

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
#import xgboost as xgb

# sklearn and model persisitence is iffy
import sklearn
if version.parse(sklearn.__version__) < version.parse('1.0.0'):
	from sklearn.experimental import enable_hist_gradient_boosting
	path = pkg_resources.resource_filename('prfect', 'clf.0.24.2.pkl')
#elif version.parse(sklearn.__version__) < version.parse('1.1.0'):
#	path = pkg_resources.resource_filename('prfect', 'clf.1.0.pkl')
#elif version.parse(sklearn.__version__) < version.parse('1.1.1'):
#	path = pkg_resources.resource_filename('prfect', 'clf.1.1.0.pkl')
else:
	path = pkg_resources.resource_filename('prfect', 'clf.1.1.1.pkl')
from sklearn.ensemble import HistGradientBoostingClassifier
clf = None

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def rint(s):
    return int(s.replace('<','').replace('>',''))

def fix_pairs(tup):
	pairs = [list(item) if item != ('1',) else ['1','1'] for item in tup]
	for i in range(2):
		if '<' in pairs[i][0]:
			pairs[i][0] = rint(pairs[i][0]) // 3 * 3 + rint(pairs[i][1]) % 3 + 1
		if '>' in pairs[i][1]:
			pairs[i][1] = rint(pairs[i][1]) // 3 * 3 + rint(pairs[i][0]) % 3 + 2

	pairs = [list(map(int,item)) for item in pairs]
	# this is to fix features that have incorrect locations by using the frame of the other end
	if pairs[0][0] % 3 != (pairs[0][1]-2) % 3:
		pairs[0][1] = pairs[0][1] // 3 * 3 + ((pairs[0][0])-1) % 3
	if pairs[1][0] % 3 != (pairs[1][1]-2) % 3:
		pairs[1][0] = pairs[1][0] // 3 * 3 + ((pairs[1][1])-2) % 3
	return tuple([tuple(map(str,pair)) for pair in pairs])

def strr(x):
    if isinstance(x, float):
        return str(round(x,3))
    else:
        return str(x)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x


def alert(args, last, curr, metrics):
	#pairs = [[last.left()+1, last.right()+1], [curr.left()+1, curr.right()+1]]
	pairs = [list(last.pairs[-1]), list(curr.pairs[0])]
	if last.strand > 0:
		pairs[0][-1] = metrics['LOC']
		pairs[-1][0] = metrics['LOC'] + 1 + metrics['DIR']
	pairs = [list(map(str, lis)) for lis in pairs] 
	feature = Feature('CDS', curr.strand, pairs, args.locus)
	
	feature.tags['ribosomal_slippage'] = [None]
	feature.tags['direction'] =  [metrics['DIR']]
	feature.tags['motif'] = [args.locus.number_motif(metrics['MOTIF'])]
	feature.tags['bases'] = [metrics['SLIPSITE']]
	feature.tags['label'] = [metrics['LABEL']]
	feature.tags['locus'] = [args.locus.name()]
	#feature.tags['location'] = [metrics['LOC']]
	#feature.tags['translation'] = feature.translation()
	if 'product' in last.tags or 'product' in curr.tags:
		feature.tags['product'] = last.tags.get('product',[]) + curr.tags.get('product',[])
	if args.format == 'feature':
		feature.write(args.outfile)

header = True
def dump(args, last, curr, metrics):
	global header
	motif = args.locus.number_motif(metrics.pop('MOTIF'))
	if header:
		keys = list(metrics.keys())
		args.outfile.print('LOCUS'.ljust(len(args.locus.name()),' '))
		args.outfile.print('\t')
		args.outfile.print('\t'.join([key.replace('R3','').ljust(5,' ') for key in map(str,keys)]))
		args.outfile.print('\t')
		args.outfile.print('MOTIF'.ljust(10))
		args.outfile.print('\t')
		args.outfile.print('GENES'.ljust(10))
		args.outfile.print('\n')
		header = False
	args.outfile.print(args.locus.name())
	args.outfile.print('\t')
	args.outfile.print('\t'.join(map(strr,metrics.values())))
	args.outfile.print('\t')
	args.outfile.print(motif.ljust(10))
	args.outfile.print('\t')
	args.outfile.print('..'.join(list(last.pairs[0])))
	args.outfile.print(',')
	args.outfile.print('..'.join(list(curr.pairs[0])))
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
	parser = argparse.ArgumentParser(description='', formatter_class=argparse.RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='feature', choices=['tabular','genbank','feature'])
	parser.add_argument('-s', '--scale', type=float, default=1.0, help="parameter to scale the MFE scores by [1.0]")
	parser.add_argument('-d', '--dump', action="store_true", help="flag to dump out the various cellular property metrics that are used for prediction")
	parser.add_argument('-m', '--model', type=str, help=argparse.SUPPRESS) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)
	
	if args.model:
		clf = pickle.load(open(args.model, 'rb'))
	else:
		clf = pickle.load(open(path, 'rb'))

	has = False	
	genbank = File(args.infile)
	for locus in genbank:
		locus.init(args)
		locus.dna = locus.dna.lower()
		args.locus = locus
		locus.args = args
		last = curr = nest = _last = _curr = None
		for curr in sorted(locus.features(include='CDS')):
			if curr.nested_inside(last):
				nest = last
			elif not curr.nested_inside(nest):
				nest = None
			for _last,_curr in set( ((last,curr),(nest,curr),(curr,nest)) ):
				best = dict()
				if _last and _curr and _last.strand == _curr.strand:
					for metrics in locus.get_metrics(_last, _curr):
						if has_prf(metrics):
							has = True
							if not best or metrics['prob'] > best['prob']:
								best = metrics
						if args.dump:
							dump(args, _last, _curr, metrics)
					if best and not args.dump:
						alert(args, _last, _curr, best)

			if curr.is_joined():
				for pairs in pairwise(curr.pairs):
					best = dict()
					loc = (int(pairs[0][-1]) + int(pairs[-1][0])) / 2
					pairs = fix_pairs(pairs)
					_last = Feature(curr.type, curr.strand, [pairs[0]], locus, curr.tags)
					_curr = Feature(curr.type, curr.strand, [pairs[1]], locus, curr.tags)
					# skip CDS with locations seperated by more than 10bp
					if abs(_curr.left() - _last.right()) < 10:
						for metrics in locus.get_metrics(_last, _curr):
							#loc = (_last.right() + _curr.left() + 4 ) / 2
							# the location of the slippery site has to be within 10bp of annotation
							metrics['LABEL'] = metrics['LOC'] - loc 
							if has_prf(metrics):
								has = True
								if not best or metrics['prob'] > best['prob']:
									best = metrics
							if args.dump:
								dump(args, _last, _curr, metrics)
						if best and not args.dump:
							alert(args, _last, _curr, best)
						elif not args.dump:
							curr.write(args.outfile)
				curr = _curr
			last = curr
if not has and not args.dump:
	print('no PRFs detected')

