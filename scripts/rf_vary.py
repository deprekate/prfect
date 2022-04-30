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

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
#from matplotlib import pyplot as plt

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
	parser.add_argument('-c', '--cluster', help='', type=int, default=0)
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
	oe = preprocessing.OrdinalEncoder()
	he = preprocessing.OneHotEncoder()

	df = pd.read_csv(args.infile, sep='\t')
	# this removes duplicate positive TYPEs that are more than 10 bases away from annotated shift
	df.loc[abs(df.I - df.CURRLEFT) > 10 ,'TYPE'] = 0

	#df['BEST'] = 0
	#df.loc[df.loc[df.TYPE==1].groupby('NAME')['PROB'].idxmin(),'BEST'] = 1
	#print(df) ; exit()


	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby('NAME')['TYPE'].any().to_frame('HAS')
	df = df.merge(has, left_on='NAME', right_index=True)
	df = df.loc[df.HAS,:]

	#df = df[df.MOTIF != "is_twoonethree"]

	#print(df.loc[df.NAME=='Abigail',["NAME",'I','TYPE','BEST']]) ; exit()

	# label encoder is prob not the best to use
	#df.loc[:,'MOTIF'] = le.fit_transform(df['MOTIF'])
	oh = pd.get_dummies(df['MOTIF'])

	take = [False] * len(df.columns)
	#for i in [1,10,11,12,16,17,18,19,20]:
	for i in [1,10,11,12,18]:
		take[i] = True

	#print(df.loc[:,take]) ; exit()


	#cluster = clusters[args.cluster-1]
	column = 'mash_k16s400c90'
	#cluster = df['CLUSTER'].unique()[args.cluster-1]
	cluster = df[column].unique()[args.cluster-1]
	colnames = df.columns
	#for model,param in [('CC','parameters_CC09.txt'), ('DP','parameters_DP09.txt')]: #'cluster in clusters:
	for model,param in [('DP','parameters_DP09.txt')]: #'cluster in clusters:
		X = df.loc[(df['MODEL']==model) & (df['PARAM']==param),:].join(oh)
		Y = df.loc[(df['MODEL']==model) & (df['PARAM']==param),['TYPE']]
		take.extend( [True] * len(oh.columns))
		for l in range(26,286,2):
			for h in range(27,286,2):
				TN = FP = FN = TP = 0
				take[l] = True
				take[h] = True
				#print(cluster, flush=True)
				#cluster = ['AikoCarson','Amok']
				#print(df.loc[:,take]) ; exit()
				X_train = X.loc[df[column] != cluster, take]
				X_test  = X.loc[df[column] == cluster, take]
				#print(X_train)  ; exit()
				Y_train = Y.loc[df[column] != cluster, :]
				Y_test  = Y.loc[df[column] == cluster, :]
		
				if X_train.empty or X_test.empty:
					continue
				clf = RandomForestClassifier(max_depth=100, random_state=0)
				clf.fit(X_train, Y_train.values.ravel())
				preds = clf.predict(X_test)
				#print(cluster)
				#print(confusion_matrix(Y_test, preds, labels=[0,1]))
				tn, fp, fn, tp = confusion_matrix(Y_test, preds, labels=[0,1]).ravel()
				TN += tn ; FP += fp ; FN += fn ; TP += tp
				#print(colored(tn, 'green'),colored(fp, 'red'),colored(fn, 'red'),colored(tp, 'green') )
		
				#print(colored(TN, 'green'),colored(FP, 'red'),colored(FN, 'red'),colored(TP, 'green') )
				precis = TP / (TP+FP) if (TP+FP) else 0
				recall = TP / (TP+FN) if (TP+FN) else 0
				f1 = (2 * precis * recall) / (precis+recall) if (precis+recall) else 0
				print(cluster, model, param, colnames[l], colnames[h], TN, FP, FN, TP, precis, recall, f1, sep='\t', flush=True)
				take[l] = False
				take[h] = False
	




