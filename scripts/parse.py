#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import sys
import re
import os


hits = dict()
with open("all.txt") as fp:
	for line in fp:
		col = line.rstrip().split('\t')
		col[3],col[4] = map(int, (col[3],col[4]))
		if col[4] / col[3] < 0.9:
			continue
		seq = col[14].upper()
		if (col[0],col[3]) in hits:
			hits[(col[0],col[3])][seq] = col
		else:
			hits[(col[0],col[3])] = {seq:col}


#print(hits['NC_004745',428])


filepath = sys.argv[1]
for filename in os.listdir(filepath):
	with open(os.path.join(filepath, filename)) as fp:
		for line in fp:
			col = line.rstrip().split('\t')
			name = eval(col[0])
			pairs = eval(col[2].replace('<','').replace('>','').replace("'", "") )
			if len(pairs) == 2:
				distance = pairs[1][0] - pairs[0][1] - 1
				if abs(distance) < 5:
					if int(col[3]) > 0:
						i = pairs[0][1] - pairs[0][0]
					else:
						i = pairs[1][1] - pairs[1][0]
					left,right = max(0,i-58),i+58
					seq = col[5].upper()
					#print(">seq" + "_" + name + "_" + str(len(seq)))
					print(col[5][left:right], flush=True, sep='\t')
					#print(name, distance, seq, flush=True, sep='\t')
					#print(seq, flush=True, sep='\t')
					for hit,value in hits[(name,len(seq))].items():
						qstart = int(value[7])
						offset = value[13][0:left-qstart].count('-')
						#print(hit[left-qstart+1+offset:right-qstart+1+offset], offset)
						a = value[13][left-qstart+1+offset:right-qstart+1+offset]
						b = value[14][left-qstart+1+offset:right-qstart+1+offset]
						#print(a)
						#print(b)
						pad = ""
						for c1,c2 in zip(a,b):
							if c1 != '-':
								print(c2, end='')
							else:
								pad += "-"
						print(pad)
					print()
					exit()
			#	print("%s.gbk -s %s..%s -d %s" % (name, pairs[0][1]-10, pairs[1][0]+10, col[3]) )
