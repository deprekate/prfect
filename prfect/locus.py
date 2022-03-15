import os
import sys
from itertools import zip_longest, chain, tee, islice
from termcolor import colored
from math import log10, exp, sqrt

from genbank.locus import Locus
from prfect.feature import Feature
import score_rbs

import numpy as np

import LinearFold as lf
from hotknots import hotknots as hk
# initialize everything first
params = os.path.dirname(hk.__file__)
hk.initialize( 'CC', os.path.join(params,"parameters_CC06.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )

def prodigal_score_rbs(seq):
	s = seq[::-1]
	score = 0
	if 'ggagga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 27
	elif 'ggagga' in (s[3:9],s[4:10]):
		score = 26
	elif 'ggagga' in (s[11:17],s[12:18]):
		score = 25
	elif 'ggagg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 24
	elif 'ggagg' in (s[3:8],s[4:9]):
		score = 23
	elif 'gagga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 22
	elif 'gagga' in (s[3:8],s[4:9]):
		score = 21
	elif 'gagga' in (s[11:16],s[12:17]) or 'ggagg' in (s[11:16],s[12:17]):
		score = 20
	elif 'ggacga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggatga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggaaga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggcgga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggggga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggtgga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggaaga' in (s[3:9],s[4:10]) or 'ggatga' in (s[3:9],s[4:10]) or 'ggacga' in (s[3:9],s[4:10]):
		score = 18
	elif 'ggtgga' in (s[3:9],s[4:10]) or 'ggggga' in (s[3:9],s[4:10]) or 'ggcgga' in (s[3:9],s[4:10]):
		score = 18
	elif 'ggaaga' in (s[11:17],s[12:18]) or 'ggatga' in (s[11:17],s[12:18]) or 'ggacga' in (s[11:17],s[12:18]):
		score = 17
	elif 'ggtgga' in (s[11:17],s[12:18]) or 'ggggga' in (s[11:17],s[12:18]) or 'ggcgga' in (s[11:17],s[12:18]):
		score = 17
	elif 'ggag' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 16
	elif 'gagg' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 16
	elif 'agga' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 15
	elif 'ggtgg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'ggggg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'ggcgg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'agg' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'gag' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'gga' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'agga' in (s[11:15],s[12:16]) or 'gagg' in (s[11:15],s[12:16]) or 'ggag' in (s[11:15],s[12:16]):
		score = 12
	elif 'agga' in (s[3:7],s[4:8]) or 'gagg' in (s[3:7],s[4:8]) or 'ggag' in (s[3:7],s[4:8]):
		score = 11
	elif 'gagga' in (s[13:18],s[14:19],s[15:20]) or 'ggagg' in (s[13:18],s[14:19],s[15:20]) or 'ggagga' in (s[13:19],s[14:20],s[15:21]):
		score = 10
	elif 'gaaga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'gatga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'gacga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'ggtgg' in (s[3:8],s[4:9]) or 'ggggg' in (s[3:8],s[4:9]) or 'ggcgg' in (s[3:8],s[4:9]):
		score = 8
	elif 'ggtgg' in (s[11:16],s[12:17]) or 'ggggg' in (s[11:16],s[12:17]) or 'ggcgg' in (s[11:16],s[12:17]):
		score = 7
	elif 'agg' in (s[11:14],s[12:15]) or 'gag' in (s[11:14],s[12:15]) or 'gga' in (s[11:14],s[12:15]):
		score = 6
	elif 'gaaga' in (s[3:8],s[4:9]) or 'gatga' in (s[3:8],s[4:9]) or 'gacga' in (s[3:8],s[4:9]):
		score = 5
	elif 'gaaga' in (s[11:16],s[12:17]) or 'gatga' in (s[11:16],s[12:17]) or 'gacga' in (s[11:16],s[12:17]):
		score = 4
	elif 'agga' in (s[13:17],s[14:18],s[15:19]) or 'gagg' in (s[13:17],s[14:18],s[15:19]) or 'ggag' in (s[13:17],s[14:18],s[15:19]):
		score = 3
	elif 'agg' in (s[13:16],s[14:17],s[15:18]) or 'gag' in (s[13:16],s[14:17],s[15:18]) or 'gga' in (s[13:16],s[14:17],s[15:18]):
		score = 2
	elif 'ggaaga' in (s[13:19],s[14:20],s[15:21]) or 'ggatga' in (s[13:19],s[14:20],s[15:21]) or 'ggacga' in (s[13:19],s[14:20],s[15:21]):
		score = 2
	elif 'ggtgg' in (s[13:18],s[14:19],s[15:20]) or 'ggggg' in (s[13:18],s[14:19],s[15:20]) or 'ggcgg' in (s[13:18],s[14:19],s[15:20]):
		score = 2
	elif 'agg' in (s[3:6],s[4:7]) or 'gag' in (s[3:6],s[4:7]) or 'gga' in (s[3:6],s[4:7]):
		score = 1
	return score


stability = [
		{'A':{'T': 2,'C': 0,'G': 0,'A': 0},'C':{'T': 0,'C':-1,'G': 3,'A': 0},'G':{'T': 1,'C': 3,'G': 0,'A': 0},'T':{'T': 0,'C': 0,'G': 1,'A': 2}}
		,
		{'A':{'T': 2,'C':-1,'G':-2,'A':-2},'C':{'T': 0,'C':-1,'G': 3,'A':-1},'G':{'T': 0,'C': 3,'G':-2,'A':-2},'T':{'T': 1,'C': 0,'G': 0,'A': 2}}
		,
		{'A':{'T': 2,'C': 0,'G': 0,'A': 0},'C':{'T': 0,'C':-1,'G': 3,'A': 0},'G':{'T': 1,'C': 3,'G': 0,'A': 0},'T':{'T': 0,'C': 0,'G': 1,'A': 2}}
		]

deltaG ={
		't':{'t':{'t':1.92,'c':2.42,'g':2.42,'a':1.25},
			 'c':{'t':3.46,'c':3.94,'g':4.30,'a':3.46},
			 'g':{'t':3.46,'c':4.30,'g':3.94,'a':3.46},
			 'a':{'t':0.04,'c':1.13,'g':1.13,'a':0.61}},
		'c':{'t':{'t':2.95,'c':3.46,'g':3.46,'a':2.30},
			 'c':{'t':4.43,'c':4.91,'g':5.27,'a':4.43},
			 'g':{'t':5.14,'c':5.98,'g':5.63,'a':5.14},
			 'a':{'t':2.30,'c':3.46,'g':3.46,'a':2.95}},
		'g':{'t':{'t':2.95,'c':3.46,'g':3.46,'a':2.30},
			 'c':{'t':5.14,'c':5.63,'g':5.98,'a':5.14},
			 'g':{'t':4.43,'c':5.27,'g':4.91,'a':4.43},
			 'a':{'t':2.30,'c':3.46,'g':3.46,'a':2.95}},
		'a':{'t':{'t':0.61,'c':1.13,'g':1.13,'a':0.04},
			 'c':{'t':3.46,'c':3.94,'g':4.30,'a':3.46},
			 'g':{'t':3.46,'c':4.30,'g':3.94,'a':3.46},
			 'a':{'t':1.25,'c':2.42,'g':2.42,'a':1.92}}
		}

def e_score(c1,c2,c3):
	global deltaG
	return exp(-deltaG[c1][c2][c3])

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

def is_ccgaaa(seq):
	if seq == 'ccgaaa':
		return True
	return False

def is_hexa(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
		return True


def is_fourtwo(seq):
	if (seq[0] == seq[1] == seq[2] == seq[3]) and (seq[4] == seq[5]):
		return True
	return False

def is_twofour(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5]):
		return True
	return False

def is_same(seq):
	return seq == len(seq) * seq[0]

def is_five(seq):
	return is_same(seq[:5])

def is_four(seq):
	return is_same(seq[:4])

def has_motif(seq):
	if is_hexa(seq) or is_five(seq) or is_fourtwo(seq) or is_twofour(seq): # or is_ccgaaa(seq):
		return True
	else:
		return False

def has(motif, seq, k):
	for i in range(len(seq)-k+1):
		if motif(seq[i:i+k]):
			return i
	return None

rbs = score_rbs.ScoreXlationInit()
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
				for codon in feature.codons():
					if codon in codons:
						codons[codon] += 1
						wobble[self.translate.wobble[codon]] += 1
		del codons['taa'] ; del codons['tag'] ; del codons['tga']
		codons['taa'],codons['tag'],codons['tga'] = 0,0,0
		del wobble['taa'] ; del wobble['tag'] ; del wobble['tga']
		total = sum(codons.values())# / pool_size
		codons = {codon:codons[codon]/total for codon in codons}
		#for k, v in sorted(codons.items(), key=lambda item: item[1]):
		#	print(k,v)
		#exit()
		#return codons

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
						counts[location[0]-1] = codons.get(codon, 0) # - count(used, codon)
						wob = self.translate.wobble.get(codon, None)
						#wounts[location[0]-1] = wobble.get(wob,ave) - count(used,wob) if codon not in ['taa','tag','tga'] else ave
						used.popleft()
						used.append(codon)
		#return counts
		#for i,vals in enumerate(grouper(counts, 3)):
		#	print(i*3+1, vals[0], vals[1], vals[2])
		for i in range(0,len(counts)-6, 3):
			print(i+1, counts[i], counts[i+1], counts[i+2])#, wounts[i], wounts[i+1], wounts[i+2])
			#for j in range(3):
			#	print(i+j+1, self.dna[i+j:i+j+3], counts[i+j])
		exit()

	def find_rbs(self):
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					i = feature.left()
					print(">seq" + str(i))
					print(self.seq(i-20, i-1, +1))
				else:
					i = feature.right()
					print(">seq" + str(i))
					print(self.seq(i+1, i+20, -1))

	def look_for_slip(self, left, right, strand):
		global rbs
		d = (1+left-right)%3-1

		if strand > 0:
			for n, i in enumerate(range(right,left,-3)):
				r = self.seq( i-19,    i,  strand)
				a = self.seq(    i+1,  i+3,  strand)
				p = self.seq(  i-2,    i,  strand)
				e = self.seq(  i-5,  i-3,  strand)
				m = self.seq(i-2+d, i+3+d, strand)
				k = self.seq(i+1,i+80, strand)
				#print(left,right, r, e,p,a, m)
				#scoring
				dist = 0 #sqrt(3*n) #10+log10(1+(right-i)/3)
				s = rbs.score_init_rbs(r, 20)
				gc = self.gc_content(k)
				# BACKWARDS
				if d < 0 and has_motif(m):
					l = lf.fold(k)[1] / gc + dist
					h = hk.fold(k.upper(), 'CC')[1] / gc + dist
					print(d, gc, left, right, r, m, [prodigal_score_rbs(r),s], e_score(*e), a, self.codon_rarity(a), l, h,  sep='\t')
					print(k)
				# FORWARD
				elif d > 0 and self.codon_rarity(a) < 0.01 and (is_four(p+a) or e==a):
					l = lf.fold(k)[1] / gc + dist
					h = hk.fold(k.upper(), 'CC')[1] / gc + dist
					print(d, gc, left, right, r, m, [prodigal_score_rbs(r),s], e_score(*e), a, self.codon_rarity(a), l, h,  sep='\t')
		return False
				#	print(colored(mfe0[1], 'red'), end='\t')

	def check_genes(self):
		#for r in self.rares:
		#	print(r, self.rares[r])
		#exit()
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
					if self.look_for_slip(left, right, _next.strand):
						sys.stderr.write(colored("ribo frameshift detected!\n", 'red') )
						print()
						print("     CDS   join(%s..%s,%s..%s)" % sum(_next.pairs, () ) )
						print()

			else:
				if _last.strand == _curr.strand and _last.frame('right') != _curr.frame('left'):
					# the left and right correspond to the min/max allowable between stop codons
					left = self.last(_curr.left()-1, _last.strand, self.stops)
					left = left + 3 if left else _curr.frame('left') - 1
					right = self.next(_last.right()-3, _last.strand, self.stops)
					right = right if right else self.length()
					if self.look_for_slip(left, right, _last.strand):
						sys.stderr.write(colored("ribo frameshift detected!\n", 'red') )
						print()
						print("     CDS   join(%s..%s,%s..%s)" % (_last.left(),_last.right(),_curr.left(),_curr.right() ) )
						print()
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



