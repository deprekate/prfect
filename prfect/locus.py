import os
import sys
from itertools import zip_longest, chain, tee, islice
from termcolor import colored
from math import log10, exp, sqrt
import pickle
import pkgutil
import pkg_resources

from genbank.locus import Locus
from prfect.feature import Feature
from prfect.functions import *
import score_rbs

import numpy as np
import pandas as pd

import LinearFold as lf
from hotknots import hotknots as hk
# initialize everything first
path = os.path.dirname(hk.__file__)
model = "CC"
param = "parameters_CC09.txt"
hk.initialize( model, os.path.join(path, param ) , os.path.join(path,"multirnafold.conf"), os.path.join(path,"pkenergy.conf") )



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

	def init(self, *args, **kwargs):
		self._rbs = score_rbs.ScoreXlationInit()
		self.backward_motifs = [is_six, is_threethree, is_fivetwo, is_twofive, is_twofour, is_threetwotwo, is_five, is_twoonefour]
		self.forward_motifs  = [is_four, is_three]
		self.motifs = self.backward_motifs + self.forward_motifs
		self.stops = ['taa','tga','tag']

	def motif_number(self, motif):
		# sklearn requires factors to be encoded into integers
		return [motif.__name__ for motif in self.motifs].index(motif)
	
	def number_motif(self, number):
		# sklearn requires factors to be encoded into integers
		return self.motifs[number]

	def score_rbs(self, rbs):
		return self._rbs.score_init_rbs(rbs,20)[0]

	def get_slips(self, last, curr):
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
			features = self.get_features(seq, d, i, j)
			if features:
				features['STOPL'] = stopL
				features['STOPR'] = stopR
				yield features
			j = j - 3

	def get_features(self, seq, d, i, j):
		features = dict()
		r  = seq[ j-23  : j-3      ]
		e0 = seq[ j-6   : j-3    ]
		p0 = seq[ j-3   : j      ]
		a0 = seq[ j     : j+3    ]
		e1 = seq[ j-6+d : j-3+d  ]
		p1 = seq[ j-3+d : j+d    ]
		a1 = seq[ j+d   : j+3+d  ]
		# features
		features['GC']     = self.gc_content()
		features['N']     = len(seq) - j - i
		features['STOPL']  = None
		features['STOPR'] = None
		features['DIR']   = d
		features['W']     = seq[ j-7 ]
		features['E0']    = e0
		features['P0']    = p0
		features['A0']    = a0
		features['Z']     = seq[ j+3 ]
		features['A0%']   = self.codon_rarity(a0)
		features['A1%']   = self.codon_rarity(a1)
		features['RBS1']  = prodigal_score_rbs(r)
		features['RBS2']  = self.score_rbs(r)
		# THIS IS TO CATCH END CASES
		#if i <= last.left()+3:
		#	return None
		if any(base not in 'acgt' for base in e1+p1+a1):
			return None
		# BACKWARDS
		elif d < 0 and self.has_backward_motif(e1+p1+a1):
			mot,prob = self.has_backward_motif(e1+p1+a1)
			features['MOTIF'] = self.motif_number(mot)
			features['PROB'] = prob
		# FORWARD
		elif d > 0 and self.has_forward_motif(e0+p0+a0) and (self.codon_rarity(a1)/self.codon_rarity(a0) > 1):
			mot,prob = self.has_forward_motif(e0+p0+a0)
			features['MOTIF'] = self.motif_number(mot)
			features['PROB'] = prob
		else:
			return None
		# deal with ambiguous bases
		seq = ''.join([base if base in 'acgt' else 'a' for base in seq])
		# ranges
		features['MODEL'] = model
		features['PARAM'] = param.replace('parameters_','').replace('.txt','')
		window = [40, 80, 120] #30,35,40,45,60,80,100] #,120]
		offset = [3] #0, 3, 6, 9, 12] #, 15]
		for w in window:
			for o in offset:
				#s = seq[ j-o-w-3 : j-o-3   ].upper().replace('T','U')
				#features['LF_%s_%s_LEFT' % (w,o)] = lf.fold(s      )[1] / len(s) / self.gc_content(s) if s else 0
				#features['HK_%s_%s_LEFT' % (w,o)] = hk.fold(s,model)[1] / len(s) / self.gc_content(s)
				s = seq[     j+o      :     j+o+w    ].upper().replace('T','U')
				features['LF_%s_%s_RIGHT' % (w,o)] = lf.fold(s      )[1] / len(s) / self.gc_content(s) if s else 0
				features['HK_%s_%s_RIGHT' % (w,o)] = hk.fold(s,model)[1] / len(s) / self.gc_content(s) if s else 0
		return features

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

