from itertools import zip_longest, chain, tee, islice

from genbank.locus import Locus
from slippery.feature import Feature

import numpy as np
import LinearFold as lf

def previous_and_next(some_iterable):
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return zip(prevs, items, nexts)

def count(dq, item):
    return sum(elem == item for elem in dq)

def is_hexa(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
		return True

def is_twofour(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5]):
		return True
	elif (seq[0] == seq[1] == seq[2] == seq[3]) and (seq[4] == seq[5]):
		return True
	return False

def is_same(seq):
	return seq == len(seq) * seq[0]

def has(motif, seq, k):
	for i in range(len(seq)-k+1):
		if motif(seq[i:i+k]):
			return True
	return False

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

	def hold(self, spacer=10):
		sites  = dict()
		for i in range(len(self.dna)-34):
			motif = self.dna[i:i+7]
			if is_hexa(motif, +1):
				mfes = []
				for j in range(0,30,3):
					seq = self.dna[ i+j+spacer : i+j+spacer+36]
					_,mfe = lf.fold( seq )
					mfes.append(mfe)	
				if has_outlier(mfes):
					for j in range(50):
						sites[i+j] = [motif] + mfes
		return sites

	def check_genes(self):
		for _last, _curr, _next in previous_and_next(sorted(self)):
			if _last is None or (_last.type != 'CDS') or (_curr.type != 'CDS'):
				pass
			elif _last.strand != _curr.strand:
				print('diff')
				pass
			elif _last.frame('right') != _curr.frame('left'):
				print('HERE')
				die()
				seq = self.seq(_last.right()-30 , _curr.left()+32)
				self.add_feature('hexa', +1, [_last.right(),_last.right()+6])
				if has(is_hexa, seq, 6):
					self.add_feature('hexa', +1, [_last.right(),_last.right()+6])
			else:
				print(_last)
				die()



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



