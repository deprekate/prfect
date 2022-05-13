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

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
from matplotlib import pyplot as plt

param_dist = {
    "max_depth": [10, None],
    "max_features": 10,
    "min_samples_split": 10,
    "bootstrap": [True, False],
    "criterion": ["gini", "entropy"],
}


def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x


if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	#parser.add_argument('clusters', type=is_valid_file, help='clusters file')
	parser.add_argument('infile', type=is_valid_file, help='input file')
	args = parser.parse_args()

	'''
	clusters = list()
	with open(args.clusters) as fp:
		for line in fp:
			cols = line.rstrip().split()
			cols = cols[1:]
			clusters.append(cols)
	'''

	le = preprocessing.LabelEncoder()
	#oe = preprocessing.OrdinalEncoder()
	#he = preprocessing.OneHotEncoder()

	df = pd.read_csv(args.infile, sep='\t')
	df = df.drop(columns=['BLANK'])
	model,param = ('DP','parameters_DP09')
	df = df.loc[(df['MODEL']==model) & (df['PARAM']==param), : ]
	# this removes duplicate positive TYPEs that are more than 10 bases away from annotated shift
	df.loc[abs(df.I - df.CURRLEFT) > 10 ,'TYPE'] = 0


	df = df[df.MOTIF != "is_twoonethree"]
	#df = df[df.MOTIF != "is_three"]
	df = df.drop( df[(df.MOTIF == "is_hexa") & (df.DIRECTION == 1)].index )

	df = df.reset_index(drop=True)

	# This will set only the best slippy sequence as TYPE==1
	#idx = df.loc[df.TYPE==1,:].groupby(['NAME','MODEL','PARAM']).PROB.idxmin()
	#df['TYPE'] = 0 ; df.loc[idx,['TYPE']] = 1
	#print(df.loc[df.TYPE==1,:]) ; exit()


	#df['LF_RATIO'] = df['L_LF_50_0'] / df['LF_50_0']
	take = ['N','RBS1','RBS2', 'a0', 'a1', 'MOTIF', 'PROB', 'LF_55_3_LEFT', 'LF_55_0_RIGHT', 'HK_30_3_LEFT','HK_30_0_RIGHT', 'LF_30_15_LEFT','LF_30_12_RIGHT']

	#take = take + list(df.columns[24:-6])
	#df = df.replace([np.inf, -np.inf], 0) 

	#df.loc[df.A == df.P, ['a0'] ] = df.loc[df.A == df.P, ['a0'] ] / 2
	#df.loc[df.A == df.E, ['a0'] ] = df.loc[df.A == df.E, ['a0'] ] / 2
	#print(df.loc[(df.NAME=='Bromden') & (df.TYPE==1) ,take]) ; exit()

	#print(df.loc[:,take+['TYPE']]) ; exit()

	#df['TYPE'] = df['TYPE'] * df['DIRECTION']
	# label encoder is prob not the best to use
	df.loc[:,'MOTIF'] = le.fit_transform(df['MOTIF'])

	#df = df.loc[df.DIRECTION==1,:]
	# this is to drop genomes that do not have a chaperone annotated
	#has = df.groupby(['NAME','DIRECTION'])['TYPE'].any().to_frame('HAS')
	has = df.groupby(['NAME','DIRECTION'])['TYPE'].any().unstack(fill_value=False).reset_index()
	df = df.merge(has) #, left_on='NAME', right_index=True)
	df = df.loc[df.HAS,:]


	TN = FP = FN = TP = 0
	#for column in ['CLUSTER','SUBCLUSTER','mash_k16s400c90','mash_k16s400c95', 'NAME']:
	for column in ['CLUSTER']:
		for cluster in df[column].unique():
			#cluster = "ClusterF"
			#print(cluster)

			#X_train = df.loc[(df[column] != cluster) & (df.DIRECTION==direction) & (df[direction]),     take     ]
			X_train = df.loc[(df[column] != cluster), take     ]
			X_test  = df.loc[(df[column] == cluster), take     ]
			Y_train = df.loc[(df[column] != cluster), ['TYPE'] ]
			Y_test  = df.loc[(df[column] == cluster), ['TYPE'] ]
			#print(X_train.join(Y_train))  ; exit()	
			#print(X_test.join(Y_test)) ; exit()	


			tot = 0 #X_test.join(df[['NAME','TYPE']]).groupby(['TYPE'])['NAME'].unique()[1].size
			
			if X_test.empty:
				continue
			#clf = RandomForestClassifier(n_estimators=100)
			#clf.fit(X_train, Y_train.values.ravel())
			clf = HistGradientBoostingClassifier(categorical_features=X_train.columns=='MOTIF', l2_regularization=0, max_iter=100).fit(X_train, Y_train.values.ravel() )
			#clf = GradientBoostingClassifier().fit(X_train, Y_train.values.ravel() )
			preds = clf.predict(X_test)
			
			tn, fp, fn, tp = confusion_matrix(Y_test, preds, labels=[0,1]).ravel()
			print(column, cluster, tn,fp,fn,tp, sep='\t', flush=True)
			'''
			for name in tem['NAME'].unique():
				tp = tem.loc[(tem.NAME==name) & (tem.TYPE==1) & (tem.PRED==1),:].shape[0]
				tn = tem.loc[(tem.NAME==name) & (tem.TYPE==0) & (tem.PRED==0),:].shape[0]
				fp = tem.loc[(tem.NAME==name) & (tem.TYPE==0) & (tem.PRED==1),:].shape[0]
				fn = tem.loc[(tem.NAME==name) & (tem.TYPE==1) & (tem.PRED==0),:].shape[0]
				for i,row in tem.loc[(tem.NAME==name) & (tem.TYPE==0) & (tem.PRED==1),:].iterrows():
					print(row.LASTLEFT, row.LASTRIGHT, row.NAME, flush=True)
				#print(column, cluster, name, tn,fp,fn,tp, sep='\t')
			'''
			#tn, fp, fn, tp = confusion_matrix(Y_test, preds, labels=[0,1]).ravel()
			#print(cluster, colored(tn, 'green'),colored(fp, 'red'),colored(fn, 'red'),colored(tp, 'green') , tot, sep='\t', flush=True)
			#print(column,cluster, model, param, TN, FP, FN, TP, precis, recall, f1, sep='\t', flush=True)
			TN += tn ; FP += fp ; FN += fn ; TP += tp
			#exit()

		#print(colored(TN, 'green'),colored(FP, 'red'),colored(FN, 'red'),colored(TP, 'green') )
		precis = TP / (TP+FP) if (TP+FP) else 0
		recall = TP / (TP+FN) if (TP+FN) else 0
		f1 = (2 * precis * recall) / (precis+recall) if (precis+recall) else 0
		print(column, precis, recall, f1, sep='\t')
	'''
		#print(tp,tn,fp,fn)
		#true = [bool(p) for p in preds]
		#print(X_test.loc[ true ,].join(Y_test.loc[true,]) )
		#sel = [bool(p) for p in Y_test.to_numpy()]
		#print(X_test.loc[ sel ,].join(Y_test.loc[sel,]) )
		#print(df.loc[ X_test.loc[ true ,].index.tolist() ,:] )
		#exit()
	'''




