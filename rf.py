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
from matplotlib import pyplot as plt

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x


if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('clusters', type=is_valid_file, help='clusters file')
	parser.add_argument('infile', type=is_valid_file, help='input file')
	args = parser.parse_args()


	clusters = list()
	with open(args.clusters) as fp:
		for line in fp:
			cols = line.rstrip().split()
			cols = cols[1:]
			clusters.append(cols)


	le = preprocessing.LabelEncoder()
	oe = preprocessing.OrdinalEncoder()
	he = preprocessing.OneHotEncoder()

	df = pd.read_csv(args.infile, sep='\t')
	# this removes duplicate positive TYPEs that are more than 10 bases away from annotated shift
	df.loc[abs(df.I - df.CURRLEFT) > 10 ,'TYPE'] = 0

	df['BEST'] = 0
	df.loc[df.loc[df.TYPE==1].groupby('NAME')['PROB'].idxmin(),'BEST'] = 1
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
	for i in [1,10,11,12,18,19,30]:
		take[i] = True

	#print(df.loc[:,take]) ; exit()

	# INBAG FOREST
	if False:
		take = [True] * len(df.columns)
		take[22] = False
		X_train = df.loc[:,take]
		Y_train = df.loc[:,'TYPE']
		clf = RandomForestClassifier(max_depth=100, random_state=0)
		clf.fit(X_train, Y_train)
		plt.barh(X_train.columns, clf.feature_importances_)
		plt.show()
		exit()


	TN = FP = FN = TP = 0
	for cluster in clusters:
		#print(cluster)
		#cluster = ['AikoCarson','Amok']
		#print(df.loc[:,take]) ; exit()
		X = df.loc[:,take].join(oh)
		X_train = X.loc[~df['NAME'].isin(cluster),:]
		X_test  = X.loc[ df['NAME'].isin(cluster),:]
		#print(X_train)  ; exit()
		Y = df.loc[:,['BEST']]
		#Y = df.loc[:,['TYPE']]
		Y_train = Y.loc[~df['NAME'].isin(cluster),:]
		Y_test  = Y.loc[ df['NAME'].isin(cluster),:]

		if X_test.empty:
			continue
		clf = RandomForestClassifier(max_depth=100, random_state=0)
		clf.fit(X_train, Y_train.values.ravel())
		preds = clf.predict(X_test)
		#print(cluster)
		'''
		print(confusion_matrix(Y_test, preds, labels=[0,1]))
		'''
		tn, fp, fn, tp = confusion_matrix(Y_test, preds, labels=[0,1]).ravel()
		TN += tn ; FP += fp ; FN += fn ; TP += tp
		#print(colored(tn, 'green'),colored(fp, 'red'),colored(fn, 'red'),colored(tp, 'green') )

	print(colored(TN, 'green'),colored(FP, 'red'),colored(FN, 'red'),colored(TP, 'green') )
	precis = TP / (TP+FP) if (TP+FP) else 0
	recall = TP / (TP+FN) if (TP+FN) else 0
	f1 = (2 * precis * recall) / (precis+recall) if (precis+recall) else 0
	print(precis, recall, f1) ; exit()
		#print(tp,tn,fp,fn)
		#true = [bool(p) for p in preds]
		#print(X_test.loc[ true ,].join(Y_test.loc[true,]) )
		#sel = [bool(p) for p in Y_test.to_numpy()]
		#print(X_test.loc[ sel ,].join(Y_test.loc[sel,]) )
		#print(df.loc[ X_test.loc[ true ,].index.tolist() ,:] )
		#exit()
	




