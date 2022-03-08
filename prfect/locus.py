import os
from itertools import zip_longest, chain, tee, islice
from termcolor import colored
from math import log10

from genbank.locus import Locus
from prfect.feature import Feature

import numpy as np

import LinearFold as lf
from hotknots import hotknots as hk
# initialize everything first
params = os.path.dirname(hk.__file__)
hk.initialize( 'CC', os.path.join(params,"parameters_CC06.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )

def previous_and_next(some_iterable):
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return zip(prevs, items, nexts)

def count(dq, item):
    return sum(elem == item for elem in dq)

def is_rbs(seq):
	motif = 'aggagg'
	match = 0
	for a,b in zip(seq,motif):
		if a==b:
			match += 1
	if match >= 5:
		return True
	else:
		return False

def is_hexa(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
		return True

def is_five(seq):
	if (seq[0] == seq[1] == seq[2] == seq[3] == seq[4]) or (seq[1] == seq[2] == seq[3] == seq[4] == seq[5]):
		return True
	return False

def is_fourtwo(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5]):
		return True
	elif (seq[0] == seq[1] == seq[2] == seq[3]) and (seq[4] == seq[5]):
		return True
	return False

def is_same(seq):
	return seq == len(seq) * seq[0]

def has_motif(seq):
	if is_hexa(seq) or is_five(seq) or is_fourtwo(seq):
		return True
	else:
		return False

def has(motif, seq, k):
	for i in range(len(seq)-k+1):
		if motif(seq[i:i+k]):
			return i
	return None

class Locus(Locus):
	def construct_feature(self):
		'''this method allows for a Feature class to be modified through inheritance in other code '''
		return Feature

	def rare_codons(self):	
		pool_size = 1000
		codons = {a+b+c : 1 for a in 'acgt' for b in 'acgt' for c in 'acgt'}
		wobble = {a+b+c : 0 for a in 'acgt' for b in 'acgt' for c in 'acgt'}
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					for codon in feature.codons():
						if codon in codons:
							codons[codon] += 1
							wobble[self.translate.wobble[codon]] += 1
		del codons['taa'] ; del codons['tag'] ; del codons['tga']
		del wobble['taa'] ; del wobble['tag'] ; del wobble['tga']
		total = sum(codons.values()) / pool_size
		codons = {codon:codons[codon]/total for codon in codons}
		total = sum(wobble.values()) / pool_size
		wobble = {wobel:wobble[wobel]/total for wobel in wobble}
		#wobble['nnn'] = 0
		from collections import deque
		ave = np.mean(list(codons.values()))
		counts = [ave] * (len(self.dna)+5)
		wounts = [ave] * (len(self.dna)+5)
		used = deque([None] * 10)
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					for codon,location in zip(feature.codons(), feature.codon_locations()):
						counts[location[0]-1] = codons.get(codon,ave) - count(used, codon)
						wob = self.translate.wobble.get(codon, None)
						wounts[location[0]-1] = wobble.get(wob,ave) - count(used,wob) if codon not in ['taa','tag','tga'] else ave
						used.popleft()
						used.append(codon)
		#for i,vals in enumerate(grouper(counts, 3)):
		#	print(i*3+1, vals[0], vals[1], vals[2])
		for i in range(0,len(counts)-6, 3):
			print(i+1, counts[i], counts[i+1], counts[i+2], wounts[i], wounts[i+1], wounts[i+2])
			#for j in range(3):
			#	print(i+j+1, self.dna[i+j:i+j+3], counts[i+j])
		exit()

	def find_slips(self):
		self.slips = [lambda : None] * self.length()
		for i in range(self.length()-5):
			seq = self.dna[i:i+6]
			self.slips[i].motif = is_hexa(seq)
		print(self.slips)

	def look_for_slip(self, left, right, strand):
		for i in range(left,right,3):
			hexa = self.dna[i:i+6]
			knot = self.dna[i+6+1:i+6+10+70]
			if has_motif(hexa):
				dist = 10+log10((right-i)/3) if strand > 0 else 10+log10((i-left)/3)
				mfe0 = lf.fold(knot)
				mfe = hk.fold(knot.upper(), 'CC')
				print(self.gc_content(), left, right, hexa, sep='\t', end='\t')
				mfe0[1] = mfe0[1] / self.gc_content() / dist
				if mfe0[1] < -5: #(-30)
					print(colored(mfe0[1], 'red'), end='\t')
				else:
					print(mfe0[1], end='\t')
				mfe[1] = mfe[1] / self.gc_content() / dist
				if mfe[1] < -5: #(-20 un normalized)
					print( colored(mfe[1], 'red'))
				else:
					print(mfe[1])

	def check_genes(self):
		self.stops = ['taa','tga','tag']
		_last = _curr = None
		for _, _, _next in previous_and_next(sorted(self)):
			#if _last is None or (_last.type != 'CDS') or (_curr.type != 'CDS'):
			if _next is None or _next.type != 'CDS':
				# ignore any features that are not CDS
				pass
			elif _last is None:
				_last = _curr
				_curr = _next
			elif len(_next.pairs) == 2:
				if _next.left() < _next.right():
					left = self.last(int(_next.pairs[1][1])-6, _next.strand, self.stops)
					left = left + 3 if left else _next.frame('right') - 1
					right = self.next(int(_next.pairs[0][0])+5, _next.strand, self.stops)
					right = right if right else self.length()
					self.look_for_slip(left, right, _next.strand)
			else:
				if _last.strand == _curr.strand and _last.frame('right') != _curr.frame('left'):
					# the left and right correspond to the min/max allowable between stop codons
					left = self.last(_curr.left()-1, _last.strand, self.stops)
					left = left + 3 if left else _curr.frame('left') - 1
					right = self.next(_last.right()-3, _last.strand, self.stops)
					right = right-3 if right else self.length()
					self.look_for_slip(left, right, _last.strand)
					'''
					if has(is_hexa, seq, 6):
						_curr.tags['__________hexa'] = True
						#self.add_feature('hexa', +1, [[_last.right(),_last.right()+6]] )
					i = has(is_rbs, seq, 6)
					if i is not None:
						_curr.tags['__________rbs'] = True
						#self.add_feature('rbs', +1, [[_last.right()+i,_last.right()+i+6]] )
					if seq:
						mfe = lf.fold(seq)
						_curr.tags['__________mfe'] = mfe[1]
						#self.add_feature('mfe=' + str(mfe[1]), +1, [[_last.right(),_curr.left()]] )
					'''	
				_last = _curr
				if _curr.right() < _next.right():
					# this is to skip CDS that completely overlap
					_curr = _next

	def heptamers(self):
		for i in range(len(self.dna)-6):
			seq = self.dna[i:i+7]
			print(i, is_hepta(seq))
		
	def mfe(self):	
		for _last, _curr, _next in previous_and_next(sorted(self)):
			if _last is None or (_last.type != 'CDS') or (_curr.type != 'CDS'):
				pass
			elif _last.strand != _curr.strand:
				pass
			elif _last.strand > 0:
				seq = self.seq(_last.right()-30 , _curr.left()+32)
				mfe = []
				for i in range(0, len(seq), 3):
					mfe.append( lf.fold(seq[i:i+30])[1] )
				Q25,Q75  = np.percentile(mfe, [75 ,25])
				IQR = Q75 - Q25
				_last.tags['mfes'] = mfe
				_last.tags['mfe'] = any(Q75 + 2*IQR > mfe)

	def write(self, outfile):
		super(Locus, self).write(outfile)



