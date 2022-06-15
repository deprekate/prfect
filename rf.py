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

import numpy as np
import pandas as pd
from sklearn.utils.class_weight import compute_sample_weight
#from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
#from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
#from matplotlib import pyplot as plt


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
	parser.add_argument('-p', '--param', type=str, default='DP03', choices=['DP03','DP09','CC06','CC09'], help="parameter set [DP03]")
	parser.add_argument('infile', type=is_valid_file, help='input file')
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	args.model = args.param[:2]

	#le = preprocessing.LabelEncoder()
	#oe = preprocessing.OrdinalEncoder()
	#he = preprocessing.OneHotEncoder()

	df = pd.read_csv(args.infile, sep='\t')

	# is no subcluster, set it to the cluster
	df.loc[df.SUBCLUSTER=='None', 'SUBCLUSTER'] =  df.loc[df.SUBCLUSTER=='None', 'CLUSTER'] + '0'

	df = df.loc[(df['MODEL']==args.model) & (df['PARAM']==args.param), : ]
	# this removes duplicate positive LABELs that are more than 10 bases away from annotated shift
	df.loc[abs(df.STOPR - 3 - df.N - df.CURRL) > 6 ,'LABEL'] = 0

	#df['RATIO'] = df['A1%'] / df['A0%'] ; df.replace([np.inf, -np.inf], 0, inplace=True)

	# This will set only the best slippy sequence as LABEL==1
	#idx = df.loc[df.LABEL==1,:].groupby(['GENOME','MODEL','PARAM']).PROB.idxmin()
	#df['LABEL'] = 0 ; df.loc[idx,['LABEL']] = 1

	take = ['DIR', 'N', 'RBS1','RBS2', 'A0%', 'A1%', 'MOTIF', 'PROB']
	#take = take + ['LF_30_3_RIGHT','HK_30_3_RIGHT']
	#take = take + ['LF_35_3_RIGHT','HK_35_3_RIGHT']
	#take = take + ['LF_40_3_RIGHT','HK_40_3_RIGHT']
	take = take + ['LF_45_3_RIGHT', 'HK_45_3_RIGHT']
	#take = take + ['LF_50_3_RIGHT', 'HK_50_3_RIGHT']
	#take = take + ['LF_60_3_RIGHT', 'HK_60_3_RIGHT']
	#take = take + ['LF_80_3_RIGHT', 'HK_80_3_RIGHT']
	take = take + ['LF_90_3_RIGHT', 'HK_90_3_RIGHT']
	#take = take + ['LF_100_3_RIGHT', 'HK_100_3_RIGHT']
	#take = take + ['LF_120_3_RIGHT', 'HK_120_3_RIGHT']

	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby(['GENOME'])['LABEL'].any().to_frame('HAS')
	df = df.merge(has, left_on='GENOME', right_index=True)
	df = df.loc[df.HAS,:]

	df['LABEL'] = df['LABEL'] * df['DIR']
	df['WEIGHT'] = compute_sample_weight(class_weight='balanced', y=df.LABEL)

	TN = FP = FN = TP = 0
	#for column in ['CLUSTER','SUBCLUSTER','MASH90','MASH95', 'GENOME']:
	for column in ['CLUSTER']:
		for cluster in df[column].unique():
			#cluster = "DW"
			cluster = None
			X_train = df.loc[(df[column] != cluster), take     ]
			X_test  = df.loc[(df[column] == cluster), take     ]
			Y_train = df.loc[(df[column] != cluster), ['LABEL'] ]
			Z_train = df.loc[(df[column] != cluster), ['WEIGHT'] ]
			Y_test  = df.loc[(df[column] == cluster), ['LABEL'] ]

			tot = X_test.join(df[['GENOME','LABEL']]).loc[lambda d: d['LABEL']!=0, 'GENOME'].nunique() #.groupby(['LABEL'])['GENOME'].unique() #[1].size
			
			if X_test.empty and cluster:
				continue
			#weights = compute_sample_weight(class_weight='balanced', y=Y_weigh.values.ravel())
			Classifier = HistGradientBoostingClassifier
			clf = Classifier(categorical_features=[c in ['MOTIF'] for c in X_train.columns], early_stopping=False, l2_regularization=10).fit(X_train, Y_train.values.ravel(), sample_weight=Z_train.values.ravel())
			pickle.dump(clf, open('all.pkl', 'wb')) ; exit()
			preds = clf.predict(X_test)
			
			tem = X_test.join(df[['LABEL']])
			tem['PRED'] = preds
			tp1 = tem.loc[(tem.LABEL==tem.PRED) & (tem.LABEL== 1),:].shape[0]
			tp2 = tem.loc[(tem.LABEL==tem.PRED) & (tem.LABEL==-1),:].shape[0]
			fp = tem.loc[(tem.LABEL==0) & (tem.PRED!=0) & (tem.DIR==tem.PRED),:].shape[0]
			fn = tem.loc[(tem.LABEL!=0) & (tem.LABEL!=tem.PRED),:].shape[0]
			tn = tem.loc[(tem.LABEL==tem.PRED) & (tem.LABEL==0),:].shape[0]

			args.outfile.print(args.model, args.param, column, cluster, tn,fp,fn,tp1+tp2, (tp1,tp2), tot, take, sep='\t')

			TN += tn ; FP += fp ; FN += fn ; TP += tp1+tp2

	precis = TP / (TP+FP) if (TP+FP) else 0
	recall = TP / (TP+FN) if (TP+FN) else 0
	f1 = (2 * precis * recall) / (precis+recall) if (precis+recall) else 0
	print(column, precis, recall, f1, sep='\t')

