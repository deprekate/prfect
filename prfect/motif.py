import os
from math import exp

from genbank.locus import Locus
from prfect.prodigal import prodigal_score_rbs
import score_rbs
import LinearFold as lf
from hotknots import hotknots as hk
# initialize everything first
path = os.path.dirname(hk.__file__)
model = "DP"
param = "parameters_DP09.txt"
hk.initialize( model, os.path.join(path, param ) , os.path.join(path,"multirnafold.conf"), os.path.join(path,"pkenergy.conf") )


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

def is_hexa(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
		return 0.001
	return None

def is_fivetwo(seq):
	if (seq[0] == seq[1] == seq[2] == seq[3] == seq[4]) and (seq[5] == seq[6]):
		return 0.001
	return None

def is_twofive(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5] == seq[6]):
		return 0.001
	return None

def is_twofour(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5]):
		return 0.004
	return None

def is_threetwotwo(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4]) and (seq[5] == seq[6]):
		return 0.004
	return None

def is_threetwo(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4]):
		return 0.0156
	return None

def is_same(seq):
	return seq == len(seq) * seq[0]

def is_six(seq):
	if is_same(seq[:6]):
		return 0.001
	return None

def is_five(seq):
	if is_same(seq[:5]):
		return 0.004
	return None

def is_four(seq):
	if is_same(seq[:4]):
		return 0.0156
	return None

def is_three(seq):
	if is_same(seq[:3]):
		return 0.0625
	return None

def is_twoonethree(seq):
	# too common
	if (seq[0] == seq[1]) and (seq[3] == seq[4] == seq[5]):
		return 0.0156
	return None

def is_twoonetwo(seq):
	# too common
	if (seq[0] == seq[1]) and (seq[3] == seq[4]):
		return 0.0625
	return None

def is_twoonefour(seq):
	if (seq[0] == seq[1]) and (seq[3] == seq[4] == seq[5] == seq[6]):
		return 0.004
	return None

def has_backward_motif(seq):
	# these are the motifs to look for
	for motif in [is_six, is_hexa, is_fivetwo, is_twofive, is_twofour, is_threetwotwo, is_five, is_twoonefour, is_twoonethree]:
		if motif(seq):
			return (motif.__name__.ljust(15), motif(seq))
	return None

def has_forward_motif(seq):
	# these are the motifs to look for
	for motif in [is_hexa]: #, is_threetwo]:
		if motif(seq):
			return (motif.__name__.ljust(15), motif(seq))
	for motif in [is_four, is_three]:
		if motif(seq[3:7]):
			return (motif.__name__.ljust(15), motif(seq[3:7]))
	return None

def has(motif, seq):
	for i in range(0,len(seq),3):
		try:
			if motif(seq[i:]):
				return True
		except:
			pass
	return False

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
		name = self.name.ljust(10)
		self.stops = ['taa','tga','tag']
		assert _last.strand==_curr.strand
		rarity = self.codon_rarity
		d = (1+_curr.left()-_last.left())%3-1
		#m1 = self.seq(_curr.left()-5  , _last.right()+5  , +1)
		#m0 = self.seq(_curr.left()-5+d, _last.right()+5+d, +1)
		left = self.last(_curr.left()-1, _last.strand, self.stops)
		left = left + 3 if left else _curr.frame('left') - 1
		right = self.next(_last.right()-3, _last.strand, self.stops)
		right = right if right else self.length()
		strand = _curr.strand
		'''
		m = m0 #if d > 0 else m1
		print(name, d, sep='\t', end='\t')
		for mot in [is_hexa, is_fivetwo, is_twofour, is_five, is_four, is_threetwotwo, is_twoonethree]:
			print(has(mot, m), end='\t')
		print(m, end='\t')
		s = self.seq(left-1+d,right+3+d)
		for mot in [is_hexa, is_fivetwo, is_twofour, is_five, is_four, is_threetwotwo, is_twoonethree]:
			print(has(mot, s), end='\t')
		print(s)
		return
		'''
		if strand > 0:
			for n, i in enumerate(range(right,left,-3)):
				r =  self.seq( i-19  , i      , strand)
				e0 = self.seq( i-5   , i-3    , strand)
				p0 = self.seq( i-2   , i      , strand)
				a0 = self.seq( i+1   , i+3    , strand)
				e1 = self.seq( i-5+d , i-3+d  , strand)
				p1 = self.seq( i-2+d , i+d    , strand)
				a1 = self.seq( i+1+d , i+3+d  , strand)
				# ratio
				r0 = rarity(a0) ; r1 = rarity(a1)
				ratio = r1/r0 if r0 else r1 * 100
				#m0 = self.seq( i-2   , i+4    , strand)
				#m1 = self.seq( i-2+d , i+4+d  , strand)
				k =  self.seq( i+7   , i+52   , strand).upper().replace('T','U')
				K =  self.seq( i+7   , i+77   , strand).upper().replace('T','U')
				#scoring
				dist = 0 #sqrt(3*n) #10+log10(1+(right-i)/3)
				gc = self.gc_content()
				#GC = self.gc_content(K)
				s = str([prodigal_score_rbs(r), self.score_rbs(r)]).ljust(8)
				nrange = range(30,95,5)
				jrange = range(0,30,3)
				# THIS IS TO CATCH END CASES
				if i <= _last.left()+3:
					pass
				# BACKWARDS
				elif d < 0 and has_backward_motif(e1+p1+a1) and lf.fold(k)[1] / len(k) / gc < -0.1:
					m,v = has_backward_motif(e1+p1+a1)
					l = lf.fold(k)[1]        / len(k) / self.gc_content(k)
					h = hk.fold(K, model)[1] / len(K) / self.gc_content(K)
					out = [name, d, gc, left, right, _last.left(), _last.right(), _curr.left(), _curr.right(), i, n, s, e0,p0,a0, r0, r1,ratio, l, h, m,v, self.v]
					out.append(model) ; out.append(param)
					for n in nrange:
						for j in jrange:
							s =  self.seq( i+1-j-n   , i-j   , strand).upper().replace('T','U')
							l = lf.fold(s)[1]        / len(s) / self.gc_content(s)
							h = hk.fold(s, model)[1] / len(s) / self.gc_content(s)
							out.append(l)
							out.append(h)
					print("\t".join([ str(rround(item)) for item in out]))
				# FORWARD
				elif d > 0 and has_forward_motif(e0+p0+a0) and rarity(a1)/rarity(a0) > 1:
					m,v = has_forward_motif(e0+p0+a0)
					l = lf.fold(k)[1]        / len(k) / self.gc_content(k)
					h = hk.fold(K, model)[1] / len(K) / self.gc_content(K)
					out = [name, d, gc, left, right, _last.left(), _last.right(), _curr.left(), _curr.right(), i, n, s, e0,p0,a0, r0, r1, ratio, l, h, m,v, self.v]
					out.append(model) ; out.append(param)
					for n in nrange:
						for j in jrange:
							s =  self.seq( i+1-j-n   , i-j   , strand).upper().replace('T','U')
							l = lf.fold(s)[1]        / len(s) / self.gc_content(s)
							h = hk.fold(s, model)[1] / len(s) / self.gc_content(s)
							out.append(l)
							out.append(h)
					print("\t".join([ str(rround(item)) for item in out]))
					#return
		#m = self.seq(_last.right() - 10, _last.right()+10, strand)
		#out = [d, self.gc_content(), left, right, m]
		#print("\t".join([ str(rround(item)) for item in out]))
		return False
				#	print(colored(mfe0[1], 'red'), end='\t')

		






