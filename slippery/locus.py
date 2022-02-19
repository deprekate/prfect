from genbank.locus import Locus
from slippery.feature import Feature

import numpy as np

def previous_and_next(some_iterable):
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return zip(prevs, items, nexts)

def count(dq, item):
    return sum(elem == item for elem in dq)

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

