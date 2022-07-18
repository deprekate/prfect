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

	take = ['N','DIR', 'RBS1','RBS2', 'MOTIF', 'A0', 'A1']
	#take = take + ['LF30R3','HK30R3']
	#take = take + ['LF40R3','HK40R3']
	#take = take + ['LF50R3', 'HK50R3']
	#take = take + ['LF60R3', 'HK60R3']
	#take = take + ['LF70R3', 'HK70R3']
	#take = take + ['LF80R3', 'HK80R3']
	#take = take + ['LF90R3', 'HK90R3']
	#take = take + ['LF100R6', 'HK100R6']
	#take = take + ['LF110R6', 'HK110R6']
	#take = take + ['LF120R6', 'HK120R6']
	take = take + list(map(lambda s : 'LF' + s, args.param.split('_'))) 
	take = take + list(map(lambda s : 'HK' + s, args.param.split('_'))) 


	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby(['GENOME'])['LABEL'].any().to_frame('HAS')
	df = df.merge(has, left_on='GENOME', right_index=True)
	#df = df.loc[df.HAS,:]

	df['DIRLABEL'] = df['DIR'] * df['LABEL']
	df['WEIGHT'] = compute_sample_weight(class_weight='balanced', y=df.DIRLABEL)

	res = df.loc[:,['GENOME','LOC','LABEL','N', 'DIR', 'HAS','MOTIF'] ]

	TN = FP = FN = TP = 0
	for column in ['CLUSTER','SUBCLUSTER','MASH90','MASH95', 'GENOME']:
		for cluster in df[column].unique():
			#cluster = None
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

			'''
			result = permutation_importance(clf, X_train, Y_train.values.ravel(), n_repeats=20, random_state=0, n_jobs=-1)
			fig, ax = plt.subplots()
			sorted_idx = result.importances_mean.argsort()
			ax.boxplot(result.importances[sorted_idx].T, vert=False, labels=np.array(take)[sorted_idx] )
			ax.set_title("Permutation Importance of each feature: " + param)
			ax.set_ylabel("Features")
			fig.tight_layout()
			plt.show()
			exit()
			'''
			# this is to be backwards compatible
			if not hasattr(clf,'feature_names_in_'):
				clf.feature_names_in_ = take

			#pickle.dump(clf, open('clf.' + sklearn.__version__ + '.pkl', 'wb')) ; exit()

			res.loc[outrows, column.lower()]  = clf.predict(X_test)
			'''
			tem = X_test.join(df[['DIRLABEL']])
			tem['PRED'] = clf.predict(X_test)
			tp = tem.loc[(tem.DIRLABEL!=0) & (tem.DIRLABEL==tem.PRED),:].shape[0]
			fn = tem.loc[(tem.DIRLABEL!=0) & (tem.DIRLABEL!=tem.PRED),:].shape[0]
			tn = tem.loc[(tem.DIRLABEL==0) & (tem.DIRLABEL==tem.PRED),:].shape[0]
			fp = tem.loc[(tem.DIRLABEL==0) & (tem.PRED!=0) & (tem.DIR==tem.PRED),:].shape[0]
			#tot = X_test.join(df[['GENOME','LABEL']]).loc[lambda d: d['LABEL']!=0, 'GENOME'].nunique() #.groupby(['LABEL'])['GENOME'].unique() #[1].size
			tot = X_test.join(df[['GENOME']]).loc[:, 'GENOME'].nunique() #.groupby(['LABEL'])['GENOME'].unique() #[1].size

			args.outfile.print(param, column, cluster, tn,fp,fn,tp, tot, take, sep='\t')
			
			TN += tn ; FP += fp ; FN += fn ; TP += tp
			'''

	mod = args.infile[5:9]
	res.to_csv(mod + '.' + args.param + '.tsv', sep='\t', index=False, na_rep=None) 
	'''
	precis = TP / (TP+FP) if (TP+FP) else 0
	recall = TP / (TP+FN) if (TP+FN) else 0
	f1 = (2 * precis * recall) / (precis+recall) if (precis+recall) else 0
	print(column, precis, recall, f1, sep='\t')
	'''
