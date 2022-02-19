import sys
import re


with open(sys.argv[1]) as fp:
	for line in fp:
		cols = line.rstrip().split('\t')
		name = eval(cols[0])
		pairs = eval(cols[2].replace("'", "") )
		if len(pairs) == 2 and abs(pairs[0][1] - pairs[1][0]) < 10:
			print("%s.gbk -s %s..%s -d %s" % (name, pairs[0][1]-10, pairs[1][0]+10, cols[3]) )
