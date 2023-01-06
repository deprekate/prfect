from math import exp

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


			

r = {'A':'T' , 'C':'G' , 'G':'C' , 'T':'A'}

stability = [
		{'A':{'T': 2,'C': 0,'G': 0,'A': 0},'C':{'T': 0,'C':-1,'G': 3,'A': 0},'G':{'T': 1,'C': 3,'G': 0,'A': 0},'T':{'T': 0,'C': 0,'G': 1,'A': 2}}
		,
		{'A':{'T': 2,'C':-1,'G':-2,'A':-2},'C':{'T': 0,'C':-1,'G': 3,'A':-1},'G':{'T': 0,'C': 3,'G':-2,'A':-2},'T':{'T': 1,'C': 0,'G': 0,'A': 2}}
		,
		{'A':{'T': 2,'C': 0,'G': 0,'A': 0},'C':{'T': 0,'C':-1,'G': 3,'A': 0},'G':{'T': 1,'C': 3,'G': 0,'A': 0},'T':{'T': 0,'C': 0,'G': 1,'A': 2}}
		]

deltaG ={
		'T':{'T':{'T':1.92,'C':2.42,'G':2.42,'A':1.25},
			 'C':{'T':3.46,'C':3.94,'G':4.30,'A':3.46},
			 'G':{'T':3.46,'C':4.30,'G':3.94,'A':3.46},
			 'A':{'T':0.04,'C':1.13,'G':1.13,'A':0.61}},
		'C':{'T':{'T':2.95,'C':3.46,'G':3.46,'A':2.30},
			 'C':{'T':4.43,'C':4.91,'G':5.27,'A':4.43},
			 'G':{'T':5.14,'C':5.98,'G':5.63,'A':5.14},
			 'A':{'T':2.30,'C':3.46,'G':3.46,'A':2.95}},
		'G':{'T':{'T':2.95,'C':3.46,'G':3.46,'A':2.30},
			 'C':{'T':5.14,'C':5.63,'G':5.98,'A':5.14},
			 'G':{'T':4.43,'C':5.27,'G':4.91,'A':4.43},
			 'A':{'T':2.30,'C':3.46,'G':3.46,'A':2.95}},
		'A':{'T':{'T':0.61,'C':1.13,'G':1.13,'A':0.04},
			 'C':{'T':3.46,'C':3.94,'G':4.30,'A':3.46},
			 'G':{'T':3.46,'C':4.30,'G':3.94,'A':3.46},
			 'A':{'T':1.25,'C':2.42,'G':2.42,'A':1.92}}
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


composition = Aminoacids(text, text2)

'''
sites = [a+b+c+d for a in 'ACGT' for b in 'ACGT' for c in 'ACGT' for d in 'ACGT']
for site in sites:
	print(site, composition.a(site))
exit()
'''


def s_score(seq):
	score = 0
	motif = 'AGGAGG'
	for i in range(6):
		if seq[i] == motif[i]:
			score += 1
	score = 0 if score < 4 else score / 3
	return score

def e_score(seq):
	global deltaG
	offset = 6
	return exp(-deltaG[seq[offset]][seq[offset+1]][seq[offset+2]])

def p_score(seq):
	global stability
	offset = 9
	m1 = 0.63
	m2 = 0.26
	s1 = s0 = 0
	for i in [0,1,2]:
		c0 = seq[offset+i]
		c1 = seq[offset+i+1]
		s0 += stability[i][c1][r[c1]]
		s1 += stability[i][c1][r[c0]]
	if seq[offset:offset+4] in ['TTTT','AAAA','AAAT']:
		return 0.3
	else:
		return (m1*s1 - m2*s0) / 3.3
'''
sites = [a+b+c+d for a in 'ACGT' for b in 'ACGT' for c in 'ACGT' for d in 'ACGT']
for site in sites:
	print(site, p_score(site))
exit()
'''
def a_score(seq):
	global composition
	offset = 12
	return composition.a(seq[offset:offset+4])



#sequence = 'CTCAGCTAGGAGGCTAGCTACGATCGATCGATCT'
sequence = 'AGGAGG' + 'TGT' + 'TTT' + 'TAAA'
sequence = '''tattcgatattacaggtctagaaggttcgcctgctattgatgtgaatgtagtgaataacg
cgacaccctccgaaacccctgctgaataaagcggggaataactattctactatgaaagtt
gtaga'''
#sequence = "AGGGGG" + "TAT" + "CTT" + "TGAC"
sequence = sequence.upper().replace('\n','').replace(' ', '').replace('T','T')

for i in range(len(sequence)-15):
	seq = sequence[i:i+16]

	E = e_score(seq)
	P = p_score(seq)
	A = a_score(seq)
	S = s_score(seq) if E+P+A > 1 else 0
	print(i, seq, S+E+P+A, ':', S, E, P, A, sep='\t')











