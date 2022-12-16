import os
import sys
from itertools import zip_longest, chain, tee, islice
from termcolor import colored
from math import log10, exp, sqrt
import pickle
import pkgutil
import pkg_resources
from math import log

from genbank.locus import Locus
from prfect.feature import Feature
from prfect.functions import *
import score_rbs

import numpy as np
import pandas as pd

import LinearFold as lf
from hotknots import hotknots as hk


def rint(s):
	return int(s.replace('<','').replace('>',''))


def rround(item, n=4):
    try:
        return round(item, n)
    except:
        try:
            return item.decode()
        except:
            return item


class Locus(Locus, feature=Feature):
	#def __init__(self,parent, *args, **kwargs):
		#self = dict(parent)
		#self.__class__ = type(parent.__class__.__name__,(self.__class__, parent.__class__),{})
		#self.__dict__ = parent.__dict__
		#self.update(dict(parent))
		#for key,value in parent.items():
		#	self[key] = value

	def init(self, args):
		self._rbs = score_rbs.ScoreXlationInit()
		self.backward_motifs = [is_six, is_threethree, is_fivetwo, is_twofive, is_twofour, is_threetwotwo, is_five, is_twoonefour]
		self.forward_motifs  = [is_four, is_three]
		self.motifs = self.backward_motifs + self.forward_motifs
		self.stops = ['taa','tga','tag']

		# initialize everything first
		path = os.path.dirname(hk.__file__)
		param = "parameters_DP09.txt"
		self.model = 'DP'
		hk.initialize( self.model, os.path.join(path, param ) , os.path.join(path,"multirnafold.conf"), os.path.join(path,"pkenergy.conf") )

	def motif_number(self, motif):
		# sklearn requires factors to be encoded into integers
		return [motif.__name__ for motif in self.motifs].index(motif)
	
	def number_motif(self, number):
		# sklearn requires factors to be encoded into integers
		return self.motifs[number]

	def score_rbs(self, rbs):
		return self._rbs.score_init_rbs(rbs,20)[0]

	def get_metrics(self, last, curr):
		assert last.strand==curr.strand , "different strands"
	
		d = (curr.left() - last.left() + 1) % 3  - 1
		if not d: return

		# this step finds the maximum possible region between two adjacent genes
		# where a frameshift could occur: before the stop codon of the first preceding
		# gene and after the furthest upstream stop codon of the following second gene
		# dont use curr.right because of stopcodon readthrough
		if last.strand > 0:
			stopL = self.last(curr.left(), self.stops, last.strand)
			stopR = self.next(last.left(), self.stops, last.strand)
			stopL = stopL if stopL else last.left() - d
			stopR = stopR if stopR else curr.right()
			if stopR >= curr.right():
				# this has to include the d
				stopR = curr.right() - d
				#stopL = curr.left()
			elif stopL <= last.left():
				stopL = last.left()
			#else:
			#	stopL = stopL - d
			stopL = stopL + 3
			stopR = stopR + 3
		else:
			stopL = self.next(curr.right(), self.stops, last.strand)
			stopR = self.last(last.right(), self.stops, last.strand)
			stopL = stopL if stopL else last.left()
			stopR = stopR if stopR else curr.right()
			if stopR >= curr.right():
				stopR = curr.right()
				#stopL = curr.left()
			elif stopL <= last.left():
				stopL = last.left()
			#else:
			#	stopR = stopR + d
			stopL = stopL
			stopR = stopR

		# the seq() method is 0-based indexed
		overlap = self.seq(stopL, stopR, curr.strand)
		if not overlap: return
		#print(stopL, stopR, overlap, d, sep='\t')
		#assert not len(overlap) % 3, "overlap error"

		# this is to pad the ends of the above maximum possible region with flanking sequence
		# in order to look for the secondary structure within it
		# there has to be a better way to do this, it may cause errors with repeated sequence
		# or short!!!!!!
		seq = self.seq(stopL-150, stopR+150, curr.strand)
		i = seq.find(overlap)
		#j = i + len(overlap)  - 3
		j = len(self.seq(stopL-150, stopL, curr.strand)) if curr.strand > 0 else len(self.seq(stopR, stopR+150, curr.strand))
		# check for slippery sequences
		n = 1
		while j > i:
			metrics = self.metrics(seq, d, j)
			if metrics:
				metrics['N'] = 3 * n
				if curr.strand > 0:
					metrics['LOC'] = stopR - metrics['N'] 
				else:
					metrics['LOC'] = stopL + metrics['N']
				yield metrics
			j = j - 3
			n += 1

	def metrics(self, seq, d, j):
		assert d , "prf direction is zero"
		metrics = dict()
		r  = seq[ j-23  : j-3    ]
		e0 = seq[ j-6   : j-3    ]
		p0 = seq[ j-3   : j      ]
		a0 = seq[ j     : j+3    ]
		e1 = seq[ j-6+d : j-3+d  ]
		p1 = seq[ j-3+d : j+0+d  ]
		a1 = seq[ j+0+d : j+3+d  ]
		# metrics
		#metrics['GC']   = self.gc_content()
		metrics['BASES'] = e0+p0+a0
		metrics['LOC']   = None
		metrics['LABEL'] = 0
		metrics['N']     = None
		metrics['DIR']   = d
		metrics['RBS1']  = prodigal_score_rbs(r)
		metrics['RBS2']  = self.score_rbs(r)
		#metrics['W']     = seq[ j-7 ]
		#metrics['E0']    = e0
		#metrics['P0']    = p0
		#metrics['A0']    = a0
		#metrics['Z']     = seq[ j+3 ]
		# THIS IS TO CATCH END CASES
		#if i <= last.left()+3:
		#	return None
		if any(base not in 'acgt' for base in e1+p1+a1):
			return None
		# BACKWARDS
		elif d < 0 and self.has_backward_motif(e1+p1+a1):
			mot,prob = self.has_backward_motif(e1+p1+a1)
			metrics['MOTIF'] = self.motif_number(mot)
			#metrics['PROB'] = prob
		# FORWARD
		elif d > 0 and self.has_forward_motif(e0+p0+a0) and self.codon_rarity(a1) >= self.codon_rarity(a0):
			mot,prob = self.has_forward_motif(e0+p0+a0)
			metrics['MOTIF'] = self.motif_number(mot)
			#metrics['PROB'] = prob
		else:
			return None
		metrics['A0']   = self.codon_rarity(a0)
		metrics['A1']   = self.codon_rarity(a1)
		# deal with ambiguous bases
		seq = ''.join([base if base in 'acgt' else 'a' for base in seq])
		# ranges
		if self.args.param:
			window = list(map(int, [i for item in self.args.param.split('_') for i in item.split('R')][::2]))
			offset = list(map(int, [self.args.param.split('R')[-1]]))
		else:
			window = [50,100] #30,40,50,60,80,90,100,120]
			offset = [0]      #0,3,6,9,12,15]
		for w in window:
			for o in offset:
				# LEFT
				#s = seq[ j-o-w-3 : j-o-3   ].upper().replace('T','U')
				#print(s, seq[j+o:j+o+w])
				#mfe = lf.fold(s)[1] if s else 0
				#metrics['LF%sL%s' % (w,o)] = mfe / len(s) / self.gc_content(s) if len(s) and self.gc_content(s) else 0
				#metrics['HK%sL%s' % (w,o)] = hk.fold(s,model)[1] / len(s) / self.gc_content(s)

				# RIGHT
				s = seq[  j+o     : j+o+w   ].upper().replace('T','U')
				#metrics['LF%sR%s' % (w,o)] = lf.fold(s      )[1] / len(s) / self.gc_content(s) if s else 0
				#metrics['HK%sR%s' % (w,o)] = hk.fold(s, self.model)[1] / len(s) / self.gc_content(s) if s else 0
				mfe = lf.fold(s)[1]
				metrics['LF%sR%s' % (w,o)] = mfe / len(s) / self.gc_content(s) if len(s) and self.gc_content(s) else 0
				#metrics['GC%s' % w] = self.gc_content(s)
				#metrics['LF%s' % w] = mfe/len(s)
				mfe = hk.fold(s, self.model)[1]
				metrics['HK%sR%s' % (w,o)] = mfe / len(s) / self.gc_content(s) if len(s) and self.gc_content(s) else 0
				#metrics['HK%s' % w] = mfe/len(s)
		return metrics

	def has_backward_motif(self, seq):
		for motif in self.backward_motifs:
			try:
				if motif(seq):
					return (motif.__name__, motif(seq))
			except:
				pass
		return None
	
	def has_forward_motif(self, seq):
		for motif in self.forward_motifs:
			try:
				if motif(seq):
					return (motif.__name__, motif(seq))
			except:
				pass
		return None

