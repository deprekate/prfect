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
		param = "parameters_CC09.txt"
		self.model = 'CC'
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
		assert last.strand==curr.strand
		d = (1+(curr.right()-2)-last.left())%3 - 1
		# this step finds the maximum possible region between two adjacent genes
		# where a frameshift could occur: before the stop codon of the first preceding
		# gene and after the furthest upstream stop codon of the following second gene
		stopL = self.last(curr.right()-6, last.strand, self.stops)
		stopL = stopL + 1 if stopL else curr.frame('left') - 1
		stopR = self.next(last.left()+2, last.strand, self.stops)
		stopR = min(curr.right()-3, stopR) + 3 if stopR else self.length()

		overlap = self.seq(stopL, stopR, curr.strand)
		if not overlap: return

		# this is to pad the ends of the above maximum possible region with flanking sequence
		# in order to look the secondary structure within it
		seq = self.seq(stopL-150, stopR+150, curr.strand)
		i = seq.find(overlap)
		j = i + len(overlap) - 3

		# check for slippery sequences
		while j > i:
			metrics = self.metrics(seq, d, i, j)
			if metrics:
				metrics['LOC'] = stopR - metrics['N'] - 2
				yield metrics
			j = j - 3

	def metrics(self, seq, d, i, j):
		metrics = dict()
		r  = seq[ j-23  : j-3      ]
		e0 = seq[ j-6   : j-3    ]
		p0 = seq[ j-3   : j      ]
		a0 = seq[ j     : j+3    ]
		e1 = seq[ j-6+d : j-3+d  ]
		p1 = seq[ j-3+d : j+d    ]
		a1 = seq[ j+d   : j+3+d  ]
		# metrics
		metrics['LOC']   = None
		metrics['LABEL'] = 0
		#metrics['GC']   = self.gc_content()
		metrics['N']     = len(seq) - j - i
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
		elif d > 0 and self.has_forward_motif(e0+p0+a0) and (self.codon_rarity(a1)/self.codon_rarity(a0) > 1):
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
		window = [50, 100] #30,40,50,60,70,80,90,100,120]
		offset = [3] #0,3,6,9,12,15]
		for w in window:
			for o in offset:
				# LEFT
				#s = seq[ j-o-w-3 : j-o-3   ].upper().replace('T','U')
				#metrics['LF%sL%s' % (w,o)] = lf.fold(s      )[1] / len(s) / self.gc_content(s) if s else 0
				#metrics['HK%sL%s' % (w,o)] = hk.fold(s,model)[1] / len(s) / self.gc_content(s)
				# RIGHT
				s = seq[     j+o      :     j+o+w    ].upper().replace('T','U')
				metrics['LF%sR%s' % (w,o)] = lf.fold(s      )[1] / len(s) / self.gc_content(s) if s else 0
				metrics['HK%sR%s' % (w,o)] = hk.fold(s, self.model)[1] / len(s) / self.gc_content(s) if s else 0
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
		#for motif in [is_hexa]: #, is_threetwo]:
		#	if motif(seq):
		#		return (motif.__name__, motif(seq))
		for motif in self.forward_motifs:
			try:
				if motif(seq[3:7]):
					return (motif.__name__, motif(seq[3:7]))
			except:
				pass
		return None

