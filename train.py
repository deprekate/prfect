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
import xgboost as xgb
import sklearn
if version.parse(sklearn.__version__) < version.parse('1.0.0'):
	from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
#from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.metrics import confusion_matrix
from sklearn.inspection import permutation_importance


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
	parser.add_argument('-p', '--param', type=str) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	parser.add_argument('infile', type=is_valid_file, help='input file')
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	
	df = pd.read_csv(args.infile, sep='\t')

	# deal with LABELS that are -1 because they are further away than where PRF annotated
	#df['LABEL'] = df['LABEL'].replace(-1, 0)

	# if subcluster none, set it to the cluster
	#df.loc[df.SUBCLUSTER=='None', 'SUBCLUSTER'] =  df.loc[df.SUBCLUSTER=='None', 'CLUSTER'] + '0'

	take = ['DIR', 'N', 'RBS1','RBS2', 'A0', 'A1', 'MOTIF']
	#take = take + ['LF50R3', 'HK50R3']
	#take = take + ['LF100R3', 'HK100R3']
	take = take + ['LF'+item for item in args.param.split('_')]
	take = take + ['HK'+item for item in args.param.split('_')]
	l2 = 1 #float(args.param)

	# this is to find genomes that do not have a chaperone annotated
	has = df.groupby(['GENOME'])['LABEL'].any().to_frame('HAS')
	df = df.merge(has, left_on='GENOME', right_index=True)
	#df = df.loc[df.HAS, :]

	df['DIRLABEL'] = df['DIR'] * df['LABEL']
	df['WEIGHT'] = compute_sample_weight(class_weight='balanced', y=df.DIRLABEL)

	X_train = df.loc[ df.HAS,     take     ]
	Y_train = df.loc[ df.HAS, ['DIRLABEL'] ].values.ravel()
	Z_train = df.loc[ df.HAS,  ['WEIGHT']  ].values.ravel()
	#Z_train = compute_sample_weight(class_weight='balanced', y=Y_train)
	
	#Z_train = compute_sample_weight(class_weight='balanced', y=df.DIRLABEL)

	#X_train.to_csv('old', index=False, sep='\t') ; exit()
	Classifier = HistGradientBoostingClassifier
	clf = Classifier(
						categorical_features=[c in ['MOTIF'] for c in X_train.columns],
						early_stopping=False,
						l2_regularization = l2
					 ).fit(
						X_train,
						Y_train,
						sample_weight=Z_train
					)

	# this is to be backwards compatible
	if not hasattr(clf,'feature_names_in_'):
		clf.feature_names_in_ = take
			
	try:
		os.makedirs( 'DP09/pkl/' + args.param )
		pass
	except:
		pass

	pickle.dump(clf, open('DP09/pkl/' + args.param + '/all.pkl', 'wb')) ; exit()
	#pickle.dump(clf, open('DP9/pkl/50R0_100R0/all.pkl', 'wb')) ; exit()
	#pickle.dump(clf, open(args.infile.split('.')[0] + '.' + args.param + '.pkl', 'wb')) ; exit()
	#pickle.dump(clf, open('cl.' + sklearn.__version__ + '.pkl', 'wb')) ; exit()
	

