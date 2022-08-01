#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import io
import os
import sys
import gzip
import faulthandler; faulthandler.enable()
import sys
import fileinput
import argparse
from argparse import RawTextHelpFormatter
from termcolor import colored
import pickle
from packaging import version

os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6
import numpy as np
import pandas as pd

import sklearn
if version.parse(sklearn.__version__) < version.parse('1.0.0'):
	from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
#from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.metrics import confusion_matrix
from sklearn.inspection import permutation_importance

#import matplotlib.pyplot as plt

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def _print(self, *items, sep=' '):
	if len(items) == 1:
		item = items[0]
		if isinstance(item, str):
			self.write(item)
		else:
			self.write(str(item))
	else: #isinstance(item, tuple):
		first = True
		for item in items:
			if not first:
				_print(self, sep)
			_print(self, item)
			first = False
		_print(self, '\n')

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('-p', '--param', type=str) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	parser.add_argument('infile', type=is_valid_file, help='input file')
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	
	df = pd.read_csv(args.infile, sep='\t')
	# if subcluster none, set it to the cluster
	df.loc[df.SUBCLUSTER=='None', 'SUBCLUSTER'] =  df.loc[df.SUBCLUSTER=='None', 'CLUSTER'] + '0'

	#take = ['GC', 'N','DIR', 'RBS1','RBS2', 'MOTIF', 'A0', 'A1']
	take = ['N','DIR', 'RBS1','RBS2', 'MOTIF', 'A0', 'A1']
	take = take + ['LF50R3', 'HK50R3']
	take = take + ['LF100R3', 'HK100R3']

	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby(['GENOME'])['LABEL'].any().to_frame('HAS')
	df = df.merge(has, left_on='GENOME', right_index=True)

	df['DIRLABEL'] = df['DIR'] * df['LABEL']
	df['WEIGHT'] = compute_sample_weight(class_weight='balanced', y=df.DIRLABEL)

	out = df.loc[:,['GENOME','LOC','LABEL','N', 'DIR', 'HAS','MOTIF'] ]

	TN = FP = FN = TP = 0
	for column in ['CLUSTER','SUBCLUSTER','MASH90','MASH95', 'GENOME']:
		for cluster in df[column].unique():
			inrows  = (df[column] != cluster) & df.HAS
			outrows = (df[column] == cluster) #  & df.HAS
			X_train = df.loc[ inrows,     take   ]
			Y_train = df.loc[ inrows, ['DIRLABEL'] ]
			Z_train = df.loc[ inrows, ['WEIGHT'] ]
			X_test  = df.loc[outrows,	  take   ]
			Y_test  = df.loc[outrows, ['DIRLABEL'] ]

			if X_test.empty and cluster:
				continue

			Classifier = HistGradientBoostingClassifier
			clf = Classifier(categorical_features=[c in ['MOTIF'] for c in X_train.columns], early_stopping=False, l2_regularization=10).fit(X_train, Y_train.values.ravel(), sample_weight=Z_train.values.ravel())

			out.loc[outrows, column.lower()]  = clf.predict(X_test)

	out.to_csv('pred.tsv', sep='\t', index=False, na_rep=None) 
