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
from itertools import tee, chain
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
#import xgboost as xgb

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

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def rint(s):
    return int(s.replace('<','').replace('>',''))

def fix_pairs(tup):
	pairs = [list(item) if item != ('1',) else ['1','1'] for item in tup]
	if '<' in pairs[0][0]:
			pairs[0][0] = rint(pairs[0][0]) // 3 * 3 + rint(pairs[0][1]) % 3 + 1
	if '>' in pairs[0][1]:
			pairs[0][1] = rint(pairs[0][1]) // 3 * 3 + rint(pairs[0][0]) % 3 + 2
	if '<' in pairs[1][0]:
			pairs[1][0] = rint(pairs[1][0]) // 3 * 3 + rint(pairs[1][1]) % 3 + 1
	if '>' in pairs[1][1]:
			pairs[1][1] = rint(pairs[1][1]) // 3 * 3 + rint(pairs[1][0]) % 3 + 2
	pairs = [list(map(int,item)) for item in pairs]
	# this is to fix features that have incorrect locations by using the frame of the other end
	if pairs[0][0] % 3 != (pairs[0][1]-2) % 3:
		pairs[0][1] = pairs[0][1] // 3 * 3 + ((pairs[0][0])-1) % 3
	if pairs[1][0] % 3 != (pairs[1][1]-2) % 3:
		pairs[1][0] = pairs[1][0] // 3 * 3 + ((pairs[1][1])-2) % 3
	return tuple([tuple(map(str,pair)) for pair in pairs])

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
	
	feature.tags['ribosomal_slippage'] =  [metrics['DIR']]
	feature.tags['motif'] = [args.locus.number_motif(metrics['MOTIF']).__name__]
	feature.tags['bases'] = [metrics['BASES']]
	feature.tags['label'] = [metrics['LABEL']]
	feature.tags['locus'] = [args.locus.name()]
	feature.tags['location'] = [metrics['LOC']]
	#feature.tags['translation'] = feature.translation()
	if 'product' in last.tags or 'product' in curr.tags:
		feature.tags['product'] = [last.tags.get('product','') , curr.tags.get('product','')]
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
	#bst = xgb.Booster({'nthread': 4})  # init model
	#bst.load_model('0001.model')  # load data
	#dtest = xgb.DMatrix(row.loc[:,clf.feature_names_in_].values, enable_categorical=True)
	#ypred = bst.predict(dtest, iteration_range=(0, bst.best_iteration + 1))
	#print(ypred)
	#return False
	#print(clf.classes_)
	prob = clf.predict_proba(row.loc[:,clf.feature_names_in_])
	#print(metrics['LABEL'], list(map(strr,prob[0])), metrics['LOC'], sep='\t')
	#return
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
		last = curr = _last = _curr = None
		for curr in sorted(locus.features(include='CDS')):
			best = dict()
			if last and last.strand == curr.strand:
				for metrics in locus.get_metrics(last, curr):
					if args.dump:
						dump(args, last, curr, metrics)
					elif has_prf(metrics):
						if not best or metrics['prob'] > best['prob']:
							best = metrics
				if best:
					alert(args, last, curr, best)
			if curr.is_joined():
				for pairs in pairwise(curr.pairs):
					best = dict()
					pairs = fix_pairs(pairs)
					_last = Feature(curr.type, curr.strand, [pairs[0]], locus, curr.tags)
					_curr = Feature(curr.type, curr.strand, [pairs[1]], locus, curr.tags)
					# skip CDS with locations seperated by more than 10bp
					if abs(_curr.left() - _last.right()) < 10:
						for metrics in locus.get_metrics(_last, _curr):
							loc = (_last.right() + _curr.left()) // 2
							# the location of the slippery site has to be within 10bp of annotation
							if loc-10 < metrics['LOC'] < loc+10: 
								metrics['LABEL'] = 1
							else:
								metrics['LABEL'] = -1
							if args.dump:
								dump(args, _last, _curr, metrics)
							elif has_prf(metrics):
								if not best or metrics['prob'] > best['prob']:
									best = metrics
						if best:
							alert(args, _last, _curr, best)
						elif not args.dump:
							curr.write(args.outfile)
							pass
				curr = _curr
			elif not best and not args.dump:
				curr.write(args.outfile)
				pass
			last = curr

