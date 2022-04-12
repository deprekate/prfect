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

	df = pd.read_csv(args.infile, sep='\t')
	# this removes duplicate positive TYPEs that are more than 10 bases away from annotated shift
	df.loc[abs(df.I - df.CURRLEFT) > 10 ,'TYPE'] = 0

	# this is to drop genomes that do not have a chaperone annotated
	has = df.groupby('NAME')['TYPE'].any().to_frame('HAS')
	df = df.merge(has, left_on='NAME', right_index=True)
	df = df.loc[df.HAS,:]

	#print(df) ; exit()

	df.loc[:,'MOTIF'] = le.fit_transform(df['MOTIF'])

	take = [False] * len(df.columns)
	#for i in [1,10,11,12,16,17,18,19,20]:
	for i in [1,10,16,17,18,19,20]:
		take[i] = True

	for cluster in clusters:
		#cluster = ['AikoCarson','Amok']
		X_train = df.loc[~df['NAME'].isin(cluster),take]
		#print(X_train) ; exit()
		X_test  = df.loc[ df['NAME'].isin(cluster),take]
		Y_train = df.loc[~df['NAME'].isin(cluster),:].TYPE
		Y_test  = df.loc[ df['NAME'].isin(cluster),:].TYPE

		if X_test.empty:
			continue
		clf = RandomForestClassifier(max_depth=100, random_state=0)
		clf.fit(X_train, Y_train)
		preds = clf.predict(X_test)
		'''
		print(cluster)
		print(confusion_matrix(Y_test, preds, labels=[0,1]))
		print()
		'''
		tn, fp, fn, tp = confusion_matrix(Y_test, preds, labels=[0,1]).ravel()
		print(tp,tn,fp,fn)
	




