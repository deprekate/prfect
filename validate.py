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
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import f1_score

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
	parser.add_argument('-g', '--genome', type=str) #, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
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
	#take = take + ['LF'+item for item in args.param.split('_')]
	#take = take + ['HK'+item for item in args.param.split('_')]
	l2 = float(args.param)

	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby(['GENOME'])['LABEL'].any().to_frame('HAS')
	df = df.merge(has, left_on='GENOME', right_index=True)

	df['DIRLABEL'] = df['DIR'] * df['LABEL']
	df['WEIGHT'] = compute_sample_weight(class_weight='balanced', y=df.DIRLABEL)

	out = df.loc[:,['GENOME','LOC','LABEL','N', 'DIR', 'HAS','MOTIF'] ]
	out['P1'] = 0
	out['P2'] = 0
	out['P3'] = 0

	enc = OneHotEncoder(handle_unknown='ignore')
	enc.fit(df[['DIRLABEL']])

	try:
		os.makedirs( 'DP09/pkl/' + args.param)
		pass
	except:
		pass

	#clu = args.param
	#args.param = '50R3_100R3'

	means = list()
	TN = FP = FN = TP = 0
	for column in ['GENOME']: #,'SUBCLUSTER','MASH95', 'GENOME']:
		#column = 'CLUSTER' #args.genome
		for cluster in df[column].unique():
			#cluster = 'A'
			#cluster = clu
			inrows  = (df[column] != cluster) & df.HAS
			outrows = (df[column] == cluster) # & df.HAS
			X_train = df.loc[ inrows,     take   ]
			Y_train = df.loc[ inrows, ['DIRLABEL'] ].values.ravel()
			#Z_train = df.loc[ inrows, [ 'WEIGHT' ] ].values.ravel()
			Z_train = compute_sample_weight(class_weight='balanced', y=Y_train)
			#
			X_test  = df.loc[outrows,	  take   ]
			Y_test  = df.loc[outrows, ['DIRLABEL'] ].values.ravel()
			#Z_test  = df.loc[outrows, [ 'WEIGHT' ] ].values.ravel()
			if X_test.empty and cluster:
				continue
			#Z_test = compute_sample_weight(class_weight='balanced', y=Y_test)

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
			#means.append(clf.score(X_test, Y_test, Z_test))
			#print(means)
			out.loc[outrows, 'P1':'P3'] = clf.predict_proba(X_test) ; continue

			try:
				os.makedirs( 'DP09/pkl/' + args.param + '/' + column)
				pass
			except:
				pass

			#pickle.dump(clf, open(args.infile.split('.')[0] + '.' + args.param + '.' + column + '_' + cluster + '.pkl', 'wb'))
			#pickle.dump(clf, open( 'pkl/' + args.param + '/all.pkl', 'wb')) ; exit()
			pickle.dump(clf, open( 'DP09/pkl/' + args.param + '/' + column + '/' + cluster + '.pkl', 'wb'))
			exit()
			#pickle.dump(clf, open( 'GENOME/' + cluster + '.pkl', 'wb')) ; exit()

	'''
	y = enc.transform(df.loc[:, ['DIRLABEL'] ]).toarray()
	yp = out.loc[:,'P1':'P3'].to_numpy()
	y_true = np.argmax(y, axis=-1)
	y_pred = np.argmax(yp, axis=-1)
	score = f1_score(y_true, y_pred, average='macro')
	print(score)
	#exit()
	'''
	#score = roc_auc_score(y, yp, multi_class='ovr', sample_weight=df.loc['WEIGHT'] )
	#print(score)
	#print(sum(means) / len(means))

	#out.to_csv('pred.tsv', sep='\t', index=False, na_rep=None) 
	out.to_csv(column + '.'+args.param + '.tsv', sep='\t', index=False, na_rep=None) 




