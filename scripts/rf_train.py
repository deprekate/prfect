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



	df = pd.read_csv(args.infile, sep='\t')
	# this removes duplicate positive TYPEs that are more than 10 bases away from annotated shift
	df.loc[abs(df.I - df.CURRLEFT) > 10 ,'TYPE'] = 0

	# this is to drop low probabilty motifs
	#df = df.loc[df.MOTIF != 'is_six', :]
	#df = df.loc[df.MOTIF != 'is_threetwotwo', :]
	#df = df.loc[df.MOTIF != 'is_twofour', :]
	df = df.loc[df.MOTIF != 'is_twoonethree', :]


	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby('NAME')['TYPE'].any().to_frame('HAS')
	df = df.merge(has, left_on='NAME', right_index=True)
	df = df.loc[df.HAS,:]

	df = df.reset_index(drop=True)
	#df['DIRECTION'] = df['DIRECTION'].astype('category')
	#df['MOTIF'] = df['MOTIF'].astype('category')

	#print(df.loc[df['TYPE'] == 1, ['MOTIF']].value_counts() ) ; exit()
	# label encoder is prob not the best to use
	#le = preprocessing.LabelEncoder()
	#df.loc[:,'MOTIF'] = le.fit_transform(df['MOTIF'])
	#enc = preprocessing.OneHotEncoder(categories='auto')
	#oh = enc.fit_transform(df['MOTIF'])
	#enc = preprocessing.OrdinalEncoder()
	#oh = pd.get_dummies(df['MOTIF'] )
	#oe = pd.DataFrame(enc.fit_transform(df.loc[ : , ['DIRECTION','MOTIF'] ] ), columns=['motif1','motif2'])

	#print(df.iloc[0:9,0:25]) ; exit()

	take = [False] * len(df.columns)
	#for i in [1,10,11,12,18,19,20]:
	for i in [0,1,10,11,12,18,19,29,21,22]:
	#for i in [1,10,11,12,16,17,18,19,181]: #  + list(range(23,34)):
		take[i] = True

	dat = pd.get_dummies(df.loc[:,take], columns=['MOTIF'])

	# INBAG FOREST
	if False:
		#take = [True] * len(df.columns) ; take[22] = False
		X_train = df.loc[:,take].join(oh)
		Y_train = df.loc[:,'TYPE']
		clf = RandomForestClassifier(max_depth=100, random_state=0)
		clf.fit(X_train, Y_train)
		plt.barh(X_train.columns, clf.feature_importances_)
		plt.show()

		exit()

	data = {'1':dict(), '-1':dict()}
	TN = FP = FN = TP = 0
	for cluster in clusters:
		#print(cluster)
		#cluster = ['AikoCarson','Amok']
		#X = df.loc[:,take].join(oh)
		#X_train = df.loc[~df['NAME'].isin(cluster),take].reset_index(drop=True)
		X_train = dat.loc[~dat['NAME'].isin(cluster), :].drop(['NAME', 'TYPE'], axis=1)
		X_test  = dat.loc[ dat['NAME'].isin(cluster), :].drop(['NAME', 'TYPE'], axis=1)
		Y_train = dat.loc[~dat['NAME'].isin(cluster),'TYPE']
		Y_test  = dat.loc[ dat['NAME'].isin(cluster),'TYPE']
		#print(X_train)  ; exit()

		if X_test.empty:
			continue
		clf = RandomForestClassifier(max_depth=100, random_state=0)
		clf.fit(X_train, Y_train)
		preds = clf.predict(X_test)
		'''
		print(confusion_matrix(Y_test, preds, labels=[0,1]))
		print()
		'''
		tn, fp, fn, tp = confusion_matrix(Y_test, preds, labels=[0,1]).ravel()
		TN += tn ; FP += fp ; FN += fn ; TP += tp

		p = X_test.join(Y_test)
		p['PRED'] = preds
		print(p)
		exit()


	print(TN, FP, FN, TP)
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
	




