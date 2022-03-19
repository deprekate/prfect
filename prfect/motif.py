import os
from math import exp

from genbank.locus import Locus
from prfect.prodigal import prodigal_score_rbs
import score_rbs
import LinearFold as lf
from hotknots import hotknots as hk
# initialize everything first
params = os.path.dirname(hk.__file__)
#hk.initialize( 'CC', os.path.join(params,"parameters_CC06.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )
#hk.initialize( 'CC', os.path.join(params,"CG_best_parameters_ISMB2007.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )
hk.initialize( 'CC', os.path.join(params, "parameters_DP09.txt" ) , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )

def rround(item, n=4):
    try:
        return round(item, n)
    except:
        try:
            return item.decode()
        except:
            return item

def shine_dalgarno(seq):
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

def is_fivetwo(seq):
	if (seq[0] == seq[1] == seq[2] == seq[3] == seq[4]) and (seq[5] == seq[6]):
		return True
	return False

def is_twofour(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5]):
		return True
	return False

def is_threetwotwo(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4]) and (seq[5] == seq[6]):
		return True
	return False

def is_same(seq):
	return seq == len(seq) * seq[0]

def is_five(seq):
	return is_same(seq[:5])

def is_four(seq):
	return is_same(seq[:4])

def has_motif(seq):
	if is_hexa(seq) or is_fivetwo(seq):
		return 0.001
	elif is_twofour(seq) or is_threetwotwo(seq) or is_five(seq): # or is_ccgaaa(seq):
		return 0.004
	return 1

class Motif(Locus):
	def __init__(self,parent, *args, **kwargs):
		#self = dict(parent)
		self.__class__ = type(parent.__class__.__name__,(self.__class__, parent.__class__),{})
		self.__dict__ = parent.__dict__
		self.update(dict(parent))
		#for key,value in parent.items():
		#	self[key] = value
		self._rbs = score_rbs.ScoreXlationInit()
		

	def score_rbs(self, rbs):
		return self._rbs.score_init_rbs(rbs,20)[0]

	def has_slip(self, _last, _curr):
		self.stops = ['taa','tga','tag']
		assert _last.strand==_curr.strand
		rarity = self.codon_rarity
		d = (1+_curr.left()-_last.left())%3-1
		left = self.last(_curr.left()-1, _last.strand, self.stops)
		left = left + 3 if left else _curr.frame('left') - 1
		right = self.next(_last.right()-3, _last.strand, self.stops)
		right = right if right else self.length()
		strand = _curr.strand
		if strand > 0:
			for n, i in enumerate(range(right,left,-3)):
				r = self.seq( i-19,    i,  strand)
				a0 = self.seq(    i+1,  i+3,  strand)
				a1 = self.seq(    i+1+d,  i+3+d,  strand)
				p = self.seq(  i-2,    i,  strand)
				p1 = self.seq(  i-2+d,    i+d,  strand)
				e = self.seq(  i-5,  i-3,  strand)
				m = self.seq(i-2+d, i+4+d, strand)
				k = self.seq(i+7,i+52, strand).upper().replace('T','U')
				K = self.seq(i+7,i+77, strand).upper().replace('T','U')
				#print(left,right, r, e,p,a, m)
				#scoring
				dist = 0 #sqrt(3*n) #10+log10(1+(right-i)/3)
				gc = self.gc_content(k)
				GC = self.gc_content(K)
				s = str([prodigal_score_rbs(r), self.score_rbs(r)]).ljust(8)
				name = self.name.ljust(10)
				# BACKWARDS
				if d < 0 and has_motif(m) < 1:
					l = lf.fold(k)
					l = l[1]/ (len(k)-l[0].count('.')) / gc
					h = hk.fold(K, 'CC')[1] / len(K) / GC
					out = [name, d, gc, left, right, n, r, m, s, a0, rarity(a0), rarity(a1),'none', l, h, self.v]
					print("\t".join([ str(rround(item)) for item in out]))
					return
				# FORWARD
				elif d > 0:
					#print(e+p, is_five(e+p), a0, a1, rarity(a0), rarity(a1))
					if a0 in self.stops:
						pass
					elif rarity(a1)/rarity(a0) > 1 and (is_four(p+a0) or is_hexa(e+p) or is_five(e+p) ):
						l = lf.fold(k)[1] / len(k) / gc
						h = hk.fold(K, 'CC')[1] / len(K) / GC
						out = [name, d, gc, left, right, n, r, m, s, a0, rarity(a0), rarity(a1), rarity(a1)/rarity(a0), l, h, self.v]
						print("\t".join([ str(rround(item)) for item in out]))
						return
		m = self.seq(_last.right() - 10, _last.right()+10, strand)
		out = [d, self.gc_content(), left, right, m]
		#print("\t".join([ str(rround(item)) for item in out]))
		return False
				#	print(colored(mfe0[1], 'red'), end='\t')

		






