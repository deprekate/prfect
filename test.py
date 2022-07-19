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
import pkgutil
import pkg_resources

os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6
import numpy as np
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

if version.parse(sklearn.__version__) < version.parse('1.0.0'):
	from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier


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

def plot_importance(clf, X_train, Y_train):
	import matplotlib.pyplot as plt
	result = permutation_importance(clf, X_train, Y_train.values.ravel(), n_repeats=20, random_state=0, n_jobs=-1)
	fig, ax = plt.subplots()
	sorted_idx = result.importances_mean.argsort()
	ax.boxplot(result.importances[sorted_idx].T, vert=False, labels=np.array(take)[sorted_idx] )
	ax.set_title("Permutation Importance of each feature: " + param)
	ax.set_ylabel("Features")
	fig.tight_layout()
	plt.show()
	exit()

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	#parser.add_argument('-p', '--param', type=str) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	parser.add_argument('infile', type=is_valid_file, help='input file')
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	
	df = pd.read_csv(args.infile, sep='\t')

	take = ['DIR', 'N', 'RBS1','RBS2', 'A0', 'A1', 'MOTIF']
	take = take + ['LF50R3', 'HK50R3']
	take = take + ['LF100R3', 'HK100R3']

	#X_train = df.loc[ df.HAS,     take     ]
	#Y_train = df.loc[ df.HAS, ['DIRLABEL'] ]
	#Z_train = df.loc[ df.HAS,  ['WEIGHT']  ]

	df = df.loc[:,['GENOME','LOC','LABEL'] + take]
	df['PRED']  = clf.predict(df.loc[:,clf.feature_names_in_])
	df.to_csv('pred.tsv', sep='\t', index=False, na_rep=None) 

	

