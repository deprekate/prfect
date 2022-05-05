#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import sys
import re
import os

from genbank.file import File
from genbank.feature import Feature
from prfect.motif import Motif



filepath = sys.argv[1]
#for filename in os.listdir(filepath):
    #genbank = File(os.path.join(filepath, filename))
if True:
	genbank = File(filepath)
	for name,locus in genbank.items():
		motif = Motif(locus)
		_last = _curr = None
		for feature in locus:
			if feature.is_type('CDS') and feature.is_joined() and '1' not in sum(feature.pairs, ()) and len(feature.pairs)==2 and int(feature.pairs[1][0])-int(feature.pairs[0][1]) < 100:
				#print(name, feature.tags['product'], sep='\t')
				a,b,c,d = map(int, sum(feature.pairs, () ))
				if abs(c-b) < 5:
					b = b - (b-a-2)%3
					c = c + (d-c-2)%3
					a,b,c,d = map(str, [a,b,c,d])
					_last = Feature(feature.type, feature.strand, [[a,b]], locus)
					_curr = Feature(feature.type, feature.strand, [[c,d]], locus)
					motif.v = True
					motif.has_slip(_last,_curr)
					_last = None
			elif feature.is_type('CDS') and len(feature.pairs)==1:
				continue
				if _last and _last.strand==feature.strand:
					motif.v = False
					motif.has_slip(_last,feature)
				_last = feature
