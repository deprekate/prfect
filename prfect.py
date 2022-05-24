#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from subprocess import Popen, PIPE, STDOUT
from types import MethodType
from termcolor import colored


sys.path.pop(0)
from genbank.feature import Feature
#import prfect
from prfect.file import File
#from prfect.motif import Motif

#def extra(self, value=None):
#	return 'extra'
#locus.rare_codons = MethodType(rare_codons, locus)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

'''
        path = pkg_resources.resource_filename('prfect', 'all.pkl')
        #data = pkgutil.get_data(__name__, "all.pkl")
        self.clf = pickle.load(open(path, 'rb'))
'''

def alert(args, label, last, curr):
	sys.stderr.write(colored("ribo frameshift detected in " + args.infile + "\n", 'red') )
	args.outfile.print("\n")
	args.outfile.print("     CDS   join(%s..%s,%s..%s)" % (last.left(), last.right(), curr.left(), curr.right()))
	args.outfile.print("           /product=%s" % last.tags['product'] )
	args.outfile.print("           /product=%s" % curr.tags['product'] )
	args.outfile.print("\n")

def dump(args, label, gc, last, curr, features):
	args.outfile.print(label)
	args.outfile.print('\t')
	args.outfile.print(gc)
	args.outfile.print('\t')
	args.outfile.print(last.left())
	args.outfile.print('\t')
	args.outfile.print(last.right())
	args.outfile.print('\t')
	args.outfile.print(curr.left())
	args.outfile.print('\t')
	args.outfile.print(curr.right())
	args.outfile.print('\t')
	args.outfile.print('\t'.join(map(str,features)))
	args.outfile.print('\n')

def _print(self, item):
	if isinstance(item, str):
		self.write(item)
	else:
		self.write(str(item))

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-d', '--dump', action="store_true")
	args = parser.parse_args()
	args.outfile.print = _print.__get__(args.outfile)

	genbank = File(args.infile)
	for name,locus in genbank.items():
		#motif = Motif(locus)
		#motif.args = args
		locus.init()
		locus.args = args
		_last = _curr = None
		for feature in locus:
			if feature.is_type('CDS') and feature.is_joined() and '1' not in sum(feature.pairs, ()) and len(feature.pairs)==2 and int(feature.pairs[1][0])-int(feature.pairs[0][1]) < 100:

				a,b,c,d = map(int, sum(feature.pairs, () ))
				if abs(c-b) < 5:
					sys.stderr.write(colored("Genome already has a joined feature:\n", 'red') )
					feature.write(sys.stderr)
					sys.stderr.write(colored("...splitting the feature into two for testing\n\n", 'red') )
					b = b - (b-a-2)%3
					c = c + (d-c-2)%3
					a,b,c,d = map(str, [a,b,c,d])
					_last = Feature(feature.type, feature.strand, [[a,b]], locus)
					_curr = Feature(feature.type, feature.strand, [[c,d]], locus)
					for slip in locus.get_slips(_last, _curr):
						if args.dump:
							dump(args, 1, locus.gc_content(), _last, _curr, slip.values())
						else:
							if False:
								alert(args, 1, _last, _curr)
					_last = None
			elif feature.is_type('CDS') and len(feature.pairs)==1:
				if _last and _last.strand==feature.strand:
					for slip in locus.get_slips(_last, feature):
						if args.dump:
							dump(args, 0, locus.gc_content(), _last, feature, slip.values())
						else:
							if False:
								alert(args, _last, feature)
				_last = feature
	

