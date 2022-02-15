#!/usr/bin/env python3

import os
import io
import sys
import gzip
import re
import argparse
import tempfile
from collections import Counter
from argparse import RawTextHelpFormatter
from itertools import zip_longest, chain
from itertools import cycle

import slippery.codons as cd
import slipper.translate as tr

import LinearFold as lf
import numpy as np

from itertools import tee, islice

def previous_and_next(some_iterable):
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return zip(prevs, items, nexts)

def count(dq, item):
    return sum(elem == item for elem in dq)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def nint(x):
	return int(x.replace('<','').replace('>',''))

def rev_comp(dna):
	a = 'acgtrykmbvdh'
	b = 'tgcayrmkvbhd'
	tab = str.maketrans(a,b)
	return dna.translate(tab)[::-1]

def has_stop(dna, strand):
	codons = ['taa','tag','tga'] if strand > 0 else ['tta','cta','tca']
	for i in range(0, len(dna), 3):
		if dna[i:i+3] in codons:
			return True
	return False

def mask(seq1, seq2):
	out1 = out2 = ''
	for tup in zip(seq1,seq2):
		if 'X' not in tup:
			out1 += tup[0]
			out2 += tup[1]
	return out1,out2

def is_hepta(seq, strand):
	if strand > 0:
		if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
			return 1
	else:
		if (seq[1] == seq[2] == seq[3]) and (seq[4] == seq[5] == seq[6]):
			return 1
	return 0

def chunker(seq, size):
	return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def grouper(iterable, n, fillvalue=None):
	args = [iter(iterable)] * n
	return zip_longest(*args, fillvalue=fillvalue)


translate = tr.Translate()

class Locus(dict):
	def __init__(self, locus, dna):
		self.locus = locus
		self.dna = dna
		self.locations = cd.Locations(self.dna)
		
		seq = self.dna + rev_comp(self.dna)
		length = len(seq)
		self.p = dict()
		self.p['a'] = seq.count('a') / length
		self.p['c'] = seq.count('c') / length
		self.p['g'] = seq.count('g') / length
		self.p['t'] = seq.count('t') / length

		self.wobble = {'gcc':'cgg', 'gct':'cgg', 'gca':'cgt', 'gcg':'cgt',
					   'aga':'tct', 'agg':'tct', 'cga':'gct', 'cgg':'gct', 'cgt':'gcg', 'cgc':'gcg',
					   'gac':'ctg', 'gat':'ctg', 
					   'aac':'ttg', 'aat':'ttg',
					   'tgc':'acg', 'tgt':'acg',
					   'gaa':'ctt', 'gag':'ctt',
					   'caa':'gtt', 'cag':'gtt',
					   'gga':'cct', 'ggg':'cct', 'ggc':'ccg', 'ggt':'ccg',
					   'cac':'gtg', 'cat':'gtg',
					   'ata':'tat', 'atc':'tag', 'att': 'tag',
					   'tta':'aat', 'ttg':'aat', 'cta':'gta', 'ctg':'gta', 'ctt':'gtg', 'ctc':'gtg',
					   'aaa':'ttt', 'aag':'ttt',
					   'atg':'tac',
					   'ttc':'aag', 'ttt':'aag',
					   'cca':'ggt', 'ccg':'ggt', 'cct':'ggg', 'ccc':'ggg',
					   'agc':'tcg', 'agt':'tcg', 'tca':'atg', 'tcg':'atg', 'tcc':'agg', 'tct':'agg',
					   'aca':'tgt', 'acg':'tgt', 'acc':'tgg', 'act':'tgg',
					   'tgg':'acc',
					   'tac':'atg', 'tat':'atg',
					   'gta':'cat', 'gtg':'cat', 'gtc':'cag', 'gtt':'cag',
					   'tag':'atc', 'taa':'att', 'tga':'act'
					   }

		self.slippage_sites = self.find_slippage()

	def seq(self, left, right):
		return self.dna[left-1 : right]

	def length(self):
		return len(self.dna)

	def pcodon(self, codon):
		codon = codon.lower()
		return self.p[codon[0]] * self.p[codon[1]] * self.p[codon[2]]

	def pstart(self):
		return self.pcodon('atg') + self.pcodon('gtg') + self.pcodon('ttg')

	def pstop(self, seq=None):
		if seq is None:
			return self.pcodon('taa') + self.pcodon('tag') + self.pcodon('tga')
		else:
			length = len(seq)
			pa = seq.count('a') / length
			pc = seq.count('c') / length
			pg = seq.count('g') / length
			pt = seq.count('t') / length
			return pt*pa*pa + pt*pa*pg + pt*pg*pa

	def find_slippage(self, spacer=10):
		sites  = dict()
		for i in range(len(self.dna)-34):
			motif = self.dna[i:i+7]
			if is_hepta(motif, +1):
				mfes = []
				for j in range(0,30,3):
					seq = self.dna[ i+j+spacer : i+j+spacer+36]
					_,mfe = lf.fold( seq )
					mfes.append(mfe)	
				if has_outlier(mfes):
					for j in range(50):
						sites[i+j] = [motif] + mfes
		return sites
					



	def rbs(self):
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					start = feature.left()+2
					feature.tags['rbs'] = self.seq(start-30,start)
				else:
					start = feature.right()
					feature.tags['rbs'] = rev_comp(self.seq(start,start+30))
	
	def rare_codons(self):	
		pool_size = 1000
		codons = {a+b+c : 1 for a in 'acgt' for b in 'acgt' for c in 'acgt'}
		wobble = {a+b+c : 0 for a in 'acgt' for b in 'acgt' for c in 'acgt'}
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					for codon in feature.codons():
						codons[codon] += 1
						wobble[self.wobble[codon]] += 1
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
		used = deque(['nnn'] * 10)
		for feature in self:
			if feature.type == 'CDS':
				if feature.strand > 0:
					for locations in feature.codon_locations():
						codon = ''.join([self.dna[loc-1] for loc in locations])
						#counts[locations[0]] = -log10(codons[codon]) if codons[codon] else 0
						#counts[locations[0]] = 1-codons[codon]-0.9 if codons[codon] else 0
						counts[locations[0]-1] = codons.get(codon,ave) - count(used, codon)
						wounts[locations[0]-1] = wobble.get(self.wobble[codon],ave) - count(used,self.wobble[codon]) if codon not in ['taa','tag','tga'] else ave
						used.popleft()
						used.append(codon)
		#for i,vals in enumerate(grouper(counts, 3)):
		#	print(i*3+1, vals[0], vals[1], vals[2])
		for i in range(0,len(counts)-6, 3):
			#print(i+1, counts[i], counts[i+1], counts[i+2], wounts[i], wounts[i+1], wounts[i+2])
			for j in range(3):
				print(i+j+1, self.dna[i+j:i+j+3], counts[i+j])
		exit()

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



	def merge(self):	
		# set dna for features and check integrity
		_last = _curr = None
		for _, _, _next in previous_and_next(sorted(self)):
			# this just mkes sure the feature locations are in the same frame
			#for i in _curr.base_locations():
			#	_curr.dna += self.dna[ i-1 : i ]
			if _last is None or (_last.type != 'CDS') or (_curr.type != 'CDS'):
				pass
			elif _curr.strand != _last.strand:
				pass
			elif _curr.frame('left') == _last.frame('right'):
				# this merges adjacent frames
				seq = self.seq(_last.right()-30 , _curr.left()+32)
				if not has_stop(seq, _last.strand):
					del self[_last]
					del self[_curr]
					_last.tags['merged'] = _last.pairs + _curr.pairs
					_last.pairs = ( (_last.left() , _curr.right()), )
					_last.tags['seq'] = seq
					_last.tags['other'] = _last.pairs
					self[_last] = True
					_curr = _next
					continue
			elif _next is None:
				pass
			elif _last.frame('right') == _next.frame('left'):
				# this merges frames broken by an embedded gene
				seq = self.seq(_last.right()-30 , _next.left()+32)
				if not has_stop(seq, _last.strand):
					del self[_last]
					del self[_next]
					_last.tags['merged'] = _last.pairs + _curr.pairs
					_last.pairs = ((_last.left() , _next.right()),)
					_last.tags['seq'] = seq
					_last.tags['pstop'] = (1-self.pstop()) ** (len(seq)/3)
					_last.tags['pstopl'] = (1-self.pstop(rev_comp(seq))) ** (len(seq)/3)
					#_last.tags['merged'] = 'true'
					_curr.tags['embedded'] = 'true'
					self[_last] = True
					_last,_curr = _curr,_last
					continue
			elif _curr.strand > 0:
				pass
			elif _curr.strand < 0:
				if abs(_curr.stop_distance()) > 100 and _curr.nearest_stop() < _last.right():
					del self[_last]
					del self[_curr]
					_last.pairs = _last.pairs + _curr.pairs
					_last.tags['joined'] = 'true'
					self[_last] = True
					_curr = _next
					continue
			_last = _curr
			_curr = _next

	def features(self, include=None, exclude=None):
		for feature in self:
			if not include or feature.type in include:
				yield feature

	def add_feature(self, key, strand, pairs):
		"""Add a feature to the factory."""
		feature = Feature(key, strand, pairs, self)
		if feature not in self:
			self[feature] = True

	def gene_coverage(self):
		''' This calculates the protein coding gene coverage, which should be around 1 '''
		cbases = tbases = 0	
		for locus in self.values():
			dna = [False] * len(locus.dna)
			seen = dict()
			for feature in locus.features(include=['CDS']):
				for i in feature.codon_locations():
					dna[i-1] = True
			cbases += sum(dna)
			tbases += len(dna)
		return 3 * cbases / tbases

	def write(self, outfile):
		outfile.write('LOCUS       ')
		outfile.write(self.locus)
		outfile.write(str(len(self.dna)).rjust(10))
		outfile.write(' bp    DNA             UNK')
		outfile.write('\n')
		outfile.write('DEFINITION  ' + self.locus + '\n')
		outfile.write('FEATURES             Location/Qualifiers\n')
		outfile.write('     source          1..')
		outfile.write(str(len(self.dna)))
		outfile.write('\n')

		for feature in sorted(self):
			feature.write(outfile)
		outfile.write('//')
		outfile.write('\n')



			
class Feature():
	def __init__(self, type_, strand, pairs, locus):
		#super().__init__(locus.locus, locus.dna)
		self.type = type_
		self.strand = strand
		# tuplize the pairs
		self.pairs = tuple([tuple(pair) for pair in pairs])
		self.locus = locus
		self.tags = dict()
		self.dna = ''
		self.partial = False

	def frame(self, end):
		if self.type != 'CDS':
			return 0
		elif end == 'right':
			return ((self.right()%3)+1) * self.strand
		elif end == 'left':
			return ((self.left()%3)+1) * self.strand

	def hypothetical(self):
		function = self.tags['product'] if 'product' in self.tags else ''
		if 'hypot'  in function or \
		   'etical' in function or \
		   'unchar' in function or \
		   ('orf' in function and 'orfb' not in function):
			return True
		else:
			return False

	def left(self):
		return int(self.pairs[0][0])
	
	def right(self):
		return int(self.pairs[-1][-1])

	def nearest_start(self):
		if self.strand > 0:
			return self.locus.locations.nearest_start(self.left(),'+')
		else:
			return self.locus.locations.nearest_start(self.right(),'-')

	def nearest_stop(self):
		if self.strand < 0:
			return self.locus.locations.nearest_stop(self.left(),'-')
		else:
			return self.locus.locations.nearest_stop(self.right(),'+')

	def start_distance(self):
		if self.strand > 0:
			return self.left() - self.nearest_start()
		else:
			return self.nearest_start() - (self.right())

	def stop_distance(self):
		if self.strand > 0:
			return self.nearest_stop() - (self.right())
		else:
			return self.left() - self.nearest_stop()

	def __str__(self):
		"""Compute the string representation of the feature."""
		return "%s\t%s\t%s\t%s" % (
				repr(self.locus.locus),
				repr(self.type),
				repr(self.pairs),
				repr(self.tags))

	def __repr__(self):
		"""Compute the string representation of the feature."""
		return "%s(%s, %s, %s, %s)" % (
				self.__class__.__name__,
				repr(self.locus),
				repr(self.type),
				repr(self.pairs),
				repr(self.tags))
	def __hash__(self):
		return hash(self.pairs)
	#def __eq__(self, other):
	#	return self.pairs == other.pairs()

	def __lt__(self, other):
		return self.left() < other.left()

	def base_locations(self, full=False):
		if full and self.partial == 'left': 
			for i in range(-((3 - len(self.dna) % 3) % 3), 0, 1):
				yield i+1
		for left,right in self.pairs:
			#left,right = map(int, [ item.replace('<','').replace('>','') for item in pair ] )
			for i in range(left,right):
				yield i

	def codon_locations(self):
		assert self.type == 'CDS'
		for triplet in grouper(self.base_locations(full=True), 3):
			if triplet[0] >= 1:
				yield triplet

	def codons(self):
		assert self.type == 'CDS'
		if self.strand > 0:
			for locations in self.codon_locations():
				yield ''.join([self.locus.dna[loc-1] for loc in locations])
		else:
			for locations in self.codon_locations():
				yield rev_comp(''.join([self.locus.dna[loc-1] for loc in locations]))
		
	
	def translation(self):
		global translate
		aa = []
		codon = ''
		first = 0 if not self.partial else len(self.dna) % 3
		for i in range(first, len(self.dna), 3):
			codon = self.dna[ i : i+3 ]
			if self.strand > 0:
				aa.append(translate.codon(codon))
			else:
				aa.append(translate.codon(rev_comp(codon)))
		if self.strand < 0:
			aa = aa[::-1]
		if aa[-1] in '#*+':
			aa.pop()
		#aa[0] = 'M'
		return "".join(aa)

	def integrity_check(self):
		seq2 = self.translation()
		if 'translation' not in self.tags:
			return 1 - ( seq2.count('#') + seq2.count('*') + seq2.count('+') ) / len(seq2)
		else:
			seq1 = self.tags['translation']
			seq1,seq2 = mask(seq1, seq2)
			seq1,seq2 = (seq1[1:], seq2[1:])
			return max(
					fuzz.ratio(seq1, seq2),
					fuzz.ratio(seq1, seq2.replace('*', 'W'))
					) / 100


	def write(self, outfile):
		outfile.write('     ')
		outfile.write( self.type.ljust(16) )
		if not self.strand > 0:
			outfile.write('complement(')
		# the pairs
		if len(self.pairs) > 1:
			outfile.write('join(')
		pairs = []
		for left, right in self.pairs:
			left = max(1,left)
			pair = str(left) + '..' + str(right+2)
			pairs.append(pair)
		outfile.write(','.join(pairs))
		if len(self.pairs) > 1:
			outfile.write(')')
		# the pairs
		if not self.strand > 0:
			outfile.write(')')
		outfile.write('\n')
		outfile.write('                     /colour=100 100 100')
		outfile.write('\n')
		for key,value in self.tags.items():
			outfile.write('                     /')
			outfile.write(str(key))
			outfile.write('=')
			outfile.write(str(value))
			outfile.write('\n')
		
		if self.type == 'CDS':
			outfile.write('                     /nearest_start=')
			#outfile.write( str(self.nearest_start) )
			outfile.write(str( self.start_distance() ))
			outfile.write('\n')
			outfile.write('                     /nearest_stop=')
			#outfile.write( str(self.nearest_stop) )
			outfile.write(str( self.stop_distance() ))
			outfile.write('\n')


	

if __name__ == "__main__":
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-f', '--format', help='Output the features in the specified format', type=str, default='tabular', choices=['tabular','genbank','fasta', 'fna','faa', 'coverage'])
	args = parser.parse_args()

	genbank = GenbankFile(args.infile)

	if args.format == 'tabular':
		for feature in genbank.features(include=['CDS']):
			args.outfile.write(str(feature))
			args.outfile.write("\n")
	elif args.format == 'coverage':
		args.outfile.write(str(genbank.gene_coverage()))
		args.outfile.write("\n")
	elif args.format == 'faa':
		for feature in genbank.features(include=['CDS']):
			args.outfile.write(">")
			args.outfile.write(feature.locus)
			args.outfile.write("[")
			args.outfile.write(feature.location.split()[1])
			args.outfile.write("]")
			args.outfile.write("\n")
			args.outfile.write(feature.translation())
			args.outfile.write("\n")

