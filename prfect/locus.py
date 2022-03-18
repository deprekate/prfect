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

def rround(item, n=4):
    try:
        return round(item, n)
    except:
        try:
            return item.decode()
        except:
            return item

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


r = {'a':'t' , 'c':'g' , 'g':'c' , 't':'a'}
stability = [
		{'a':{'t': 2,'c': 0,'g': 0,'a': 0},'c':{'t': 0,'c':-1,'g': 3,'a': 0},'g':{'t': 1,'c': 3,'g': 0,'a': 0},'t':{'t': 0,'c': 0,'g': 1,'a': 2}}
		,
		{'a':{'t': 2,'c':-1,'g':-2,'a':-2},'c':{'t': 0,'c':-1,'g': 3,'a':-1},'g':{'t': 0,'c': 3,'g':-2,'a':-2},'t':{'t': 1,'c': 0,'g': 0,'a': 2}}
		,
		{'a':{'t': 2,'c': 0,'g': 0,'a': 0},'c':{'t': 0,'c':-1,'g': 3,'a': 0},'g':{'t': 1,'c': 3,'g': 0,'a': 0},'t':{'t': 0,'c': 0,'g': 1,'a': 2}}
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
text = '''Ala1	A	1	TGC	GCT,GCA,GCG	3250	4.55	0.2571	0.0014
Ala2	A	2	GGC	GCC	617	0.86	0.2571	0.0073
Arg2	R	3	ACG	CGT,CGC,CGA	4752	6.65	0.2568	0.0009
Arg3	R	4	CCG	CGG	639	0.89	0.2568	0.007
Arg4	R	5	TCT	AGA	867	1.21	0.2568	0.0052
Arg5	R	6	CCT	AGG	420	0.59	0.2568	0.0107
Asn	N	7	GTT	AAC,AAT	1193	1.67	0.2570	0.0038
Asp 1	D	8	GTC	GAC,GAT	2396	3.35	0.2570	0.0019
Cys	C	9	GCA	TGC,TGT	1587	2.22	0.2570	0.0028
Gln1	Q	10	TTG	CAA	764	1.07	0.2569	0.0059
Gln2	Q	11	CTG	CAG	881	1.23	0.2569	0.0051
Glu2	E	12	TTC	GAA,GAG	4717	6.60	0.2569	0.0009
Gly1	G	13	CCC	GGG	1068.5	1.49	0.2572	0.0042
Gly2	G	14	TCC	GGA,GGG	1068.5	1.49	0.2572	0.0042
Gly3	G	15	GCC	GGC,GGT	4359	6.10	0.2572	0.001
His	H	16	GTG	CAC,CAT	639	0.89	0.2569	0.007
Ile1	I	17	GAT	ATC,ATT	1737	2.43	0.2570	0.0026
Ile2	I	18	CAT	ATA	1737	2.43	0.2570	0.0026
Leu 1	L	19	CAG	CTG	4470	6.25	0.2570	0.001
Leu 2	L	20	GAG	CTC,CTT	943	1.32	0.2570	0.0048
Leu 3	L	21	TAG	CTA,CTG	666	0.93	0.2570	0.0067
Leu 4	L	22	CAA	TTG	1913	2.68	0.2570	0.0023
Leu 5	L	23	TAA	TTA,TTG	1031	1.44	0.2570	0.0043
Lys	K	24	TTT	AAA,AAG	1924	2.69	0.2569	0.0023
Met f1	M	25	CAT	ATG	1211	1.69	0.2569	0.0037
Met f2	M	26	CAT	ATG	715	1.00	0.2569	0.0063
Met m	M	27	CAT	ATG	706	0.99	0.2569	0.0064
Phe	F	28	GAA	TTC,TTT	1037	1.45	0.2568	0.0043
Pro1	P	29	CGG	CCG	900	1.26	0.2570	0.005
Pro2	P	30	GGG	CCC,CCT	720	1.01	0.2570	0.0063
Pro3	P	31	TGG	CCA,CCT,CCG	581	0.81	0.2570	0.0077
Sec	X	32	TCA	TGA	219	0.31	0.2575	0.0204
Ser1	S	33	TGA	TCA,TCT,TCG	1296	1.81	0.2571	0.0035
Ser2	S	34	CGA	TCG	344	0.48	0.2571	0.0131
Ser3	S	35	GCT	AGC,AGT	1408	1.97	0.2571	0.0032
Ser5	S	36	GGA	TCC,TCT	764	1.07	0.2571	0.0059
Thr1	T	37	GGT	ACC,ACT	104	0.15	0.2570	0.0434
Thr2	T	38	CGT	ACG	541	0.76	0.2570	0.0083
Thr3	T	39	GGT	ACC,ACT	1095	1.53	0.2570	0.0041
Thr4	T	40	TGT	ACA,ACT,ACG	916	1.28	0.2570	0.0049
Trp	W	41	CCA	TGG	943	1.32	0.2567	0.0046
Tyr1	Y	42	GTA	TAC,TAT	769	1.08	0.2568	0.0058
Tyr2	Y	43	GTA	TAC,TAT	1261	1.76	0.2568	0.0036
Val1	V	44	TAC	GTA,GTG,GTT	3840	5.37	0.2570	0.0012
Val2 A	V	45	GAC	GTC,GTT	630	0.88	0.2570	0.0072
Val2 B	V	46	GAC	GTC,GTT	635	0.89	0.2570	0.0071
RF1	X	47		TAA,TAG	1200	1.68	0.3947	0.0003
RF2	X	48		TAA,TGA	6000	8.39	0.3947	0.0001'''
text2 = '''TTT	28	22,23	GTT	44,45,46	
TTC	28	9,22,23,36,42,43	GTC	45,46	2,8,15,44
TTG	22,23	28,34,41	GTG	44	13,45,46
TTA	23	22,28,32,33	GTA	44	1,12,14,45,46
TCT	33,36	34	GCT	1	2
TCC	36	9,28,33,34,42,43	GCC	2	1,8,15,45,46
TCG	33	22,36,41	GCG	1	2,13
TCA	33	23,32,34,36	GCA	1	2,12,14,44
TGT	9	32,41	GGT	15	13,14
TGC	9	28,32,36,41,42,43	GGC	15	2,8,13,14,45,46
TGG	41	9,22,32,34	GGG	13,14	15
TGA	32,48	9,23,33,41	GGA	14	1,12,13,15,44
TAT	42,43		GAT	8	12
TAC	42,43	9,28,36	GAC	8	2,12,15,45,46
TAG	47	22,34,41,42,43	GAG	12	8,13
TAA	47,48	23,32,33,42,43	GAA	12	1,8,14,44
CTT	20	3,19,21	ATT	17	18,25,26,27
CTC	20	16,19,21,30	ATC	17	7,18,25,26,27,35,37,39
CTG	19,21	4,11,20,29	ATG	27	6,17,18,25,26,38
CTA	21	10,19,20,31	ATA	18	5,17,24,25,26,27,40
CCT	30,31	3,29	ACT	37,39,40	38
CCC	22	16,20,29,31	ACC	37,39	7,17,35,38,40
CCG	29,31	4,11,19,30	ACG	38,40	6,18,25,26,27,37,39
CCA	31	10,21,29,30	ACA	40	5,24,37,38,39
CGT	3	4	AGT	35	5,6
CGC	3	4,16,20,30	AGC	35	5,6,7,17,37,39
CGG	4	3,11,19,29	AGG	6	5,18,25,26,27,35,38
CGA	3	4,10,21,31	AGA	5	6,24,35,40
CAT	16	3,10,11	AAT	7	24
CAC	16	10,11,20,30	AAC	7	17,24,35,37,39
CAG	11	4,10,16,19,29	AAG	24	6,7,18,25,26,27,38
CAA	10	11,16,21,31	AAA	24	5,7,40'''
class Aminoacids(list):
	def __init__(self, text1, text2):
		self.codons = {}
		self.cog = {}
		self.ncog = {}
		for line in text1.split('\n'):
			cols = line.split('\t')
			cols[5] = float(cols[5])
			cols[8] = float(cols[8])
			self.append(cols)
			for codon in cols[4].split(','):
				self.codons[codon] = int(cols[2])
		for line in text2.split('\n'):
			cols = line.replace(',',' ').split('\t')
			self.cog[cols[0]] = list(map(int, cols[1].split()))
			self.ncog[cols[0]] = list(map(int, cols[2].split()))
			self.cog[cols[3]] = list(map(int, cols[4].split()))
			self.ncog[cols[3]] = list(map(int, cols[5].split()))
	def ratio(self, codon):
		ncog = cog = 0
		for num in self.ncog[codon]:
			#ncog += self[num-1][8]**-1
			ncog += self[num-1][5]
		for num in self.cog[codon]:
			#cog += self[num-1][8]**-1
			cog += self[num-1][5]
		return ncog / cog
	def cratio(self, site):
		codon1 = site[1:  ]
		codon0 = site[ :-1]
		#return self[self.codons[codon1]-1][8]**-1 / self[self.codons[codon0]-1][8]**-1
		return self[self.codons[codon1]-1][5] / self[self.codons[codon0]-1][5]
	def get(self, value, i):
		if type(value) is int:
			return self[value-1][i]
		else:
			return self[self.codons[value]-1][i]
	def molecules(self, value):
		return self.get(value,5)
	def fraction(self, value):
		return self.get(value,6)
	def diffusion(self, value):
		return self.get(value,7)
	def arrival(self, value):
		return self.get(value,8)
	def a(self, site):
		a0 = a1 = 0
		codon = site[:-1]
		a0 = self.ratio(codon) / 14.94
		a1 = self.cratio(site) / 13.8
		if codon in ['TAG','TGA']:
			a0 = a1 = 0.9
		elif codon in ['TAA']:
			a0 = a1 = 0.6
		return a0/14.94 + a1/13.8
composition = Aminoacids(text, text2)

def e_score(c1,c2,c3):
	global deltaG
	return round(exp(-deltaG[c1][c2][c3]), 4)

def p_score(p,p1):
	global stability
	m1 = 0.63
	m2 = 0.26
	s1 = s0 = 0
	for i in [0,1,2]:
		c0 = p[i]
		c1 = p1[i]
		s0 += stability[i][c1][r[c1]]
		s1 += stability[i][c1][r[c0]]
	return round((m1*s1 - m2*s0) / 3.3, 5)

def a_score(seq):
    global composition
    offset = 12
    return composition.a(seq[offset:offset+4])


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

def is_fivetwo(seq):
	if (seq[0] == seq[1] == seq[2] == seq[3] == seq[4]) and (seq[5] == seq[6]):
		return True
	return False

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
	if is_hexa(seq) or is_fivetwo(seq):
		return 0.001
	elif is_twofour(seq): # or is_ccgaaa(seq):
		return 0.004
	else:
		return None

def has(motif, seq, k):
	for i in range(len(seq)-k+1):
		if motif(seq[i:i+k]):
			return i
	return None

rbs = score_rbs.ScoreXlationInit()
class Locus(Locus, feature=Feature):

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
		for k, v in sorted(codons.items(), key=lambda item: item[1]):
			print(k,v)
		exit()
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
		rarity = self.codon_rarity
		d = (1+left-right)%3-1

		if strand > 0:
			for n, i in enumerate(range(right,left,-3)):
				r = self.seq( i-19,    i,  strand)
				a0 = self.seq(    i+1,  i+3,  strand)
				a1 = self.seq(    i+2,  i+4,  strand)
				p = self.seq(  i-2,    i,  strand)
				p1 = self.seq(  i-2+d,    i+d,  strand)
				e = self.seq(  i-5,  i-3,  strand)
				m = self.seq(i-2+d, i+4+d, strand)
				k = self.seq(i+4,i+84, strand)
				#print(left,right, r, e,p,a, m)
				#scoring
				dist = 0 #sqrt(3*n) #10+log10(1+(right-i)/3)
				s = str([prodigal_score_rbs(r), rbs.score_init_rbs(r, 20)[0]]).ljust(8)
				gc = self.gc_content(k)
				# BACKWARDS
				if d < 0 and has_motif(m):
					l = lf.fold(k)[1] / gc + dist
					h = hk.fold(k.upper(), 'CC')[1] / gc + dist
					if (l+h)/2 < -45:
						v = has_motif(m) * exp(h/80)
						out = [d, gc, left, right, r, m, v, s, e_score(*e), p_score(p,p1), a0, rarity(a0), rarity(a1),'none', l, h]
						print("\t".join([ str(rround(item)) for item in out]))
				# FORWARD
				elif d > 0:
					if a0 in self.stops:
						pass
					elif rarity(a1)/rarity(a0) > 1 and (is_four(p+a0)): # or a0=='ccc'):
						l = lf.fold(k)[1] / gc + dist
						h = hk.fold(k.upper(), 'CC')[1] / gc + dist
						if (l+h) / 2 < -45:
							v = 0.0156 * exp(h/80) * rarity(a0) / rarity(a1)
							out = [d, gc, left, right, r, m, v, s, e_score(*e), p_score(p,p1), a0, rarity(a0), rarity(a1), rarity(a1)/rarity(a0), l, h]
							print("\t".join([ str(rround(item)) for item in out]))
						#print(k)
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



