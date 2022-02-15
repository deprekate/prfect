
class Translate:
	def __init__(self):
		nucs = 'acgt'
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		amino_acids = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV#Y+YSSSS*CWCLFLF'
		self.translate = dict(zip(codons, amino_acids))
		ambi = """aarK aayN acrT acyT ackT acmT acbT acvT acdT achT agrR agyS atyI atmI athI 
		carQ cayH ccrP ccyP cckP ccmP ccbP ccvP ccdP cchP cgrR cgyR cgkR cgmR cgbR cgvR cgdR 
		cghR ctrL ctyL ctkL ctmL ctbL ctvL ctdL cthL garE gayD gcrA gcyA gckA gcmA gcbA gcvA 
		gcdA gchA ggrG ggyG ggkG ggmG ggbG ggvG ggdG gghG gtrV gtyV gtkV gtmV gtbV gtvV gtdV 
		gthV tar* tayY tcrS tcyS tckS tcmS tcbS tcvS tcdS tchS tgyC ttrL ttyF tra* ytaL ytgL 
		ytrL mgaR mggR mgrR"""
		for item in ambi.split():
			self.translate[item[0:3]] = item[-1]
		self.amino_acids = sorted(set(amino_acids))

	def rev_comp(self, seq):
		seq_dict = {'a':'t','t':'a','g':'c','c':'g',
					'n':'n',
					'r':'y','y':'r','s':'s','w':'w','k':'m','m':'k',
					'b':'v','v':'b','d':'h','h':'d'}
		return "".join([seq_dict[base] for base in reversed(seq)])

	def codon(self, codon):
		if len(codon) == 3:
			return self.translate.get(codon.lower(), 'X')
		else:
			return ''

	def sequence(self, seq, strand):
		aa = ''
		if strand > 0:
			for i in range(0, len(seq), 3):
				aa += self.codon(seq[i:i+3])
			return aa
		else:
			for i in range(0, len(seq), 3):
				aa += self.codon(self.rev_comp(seq[i:i+3]))
			return aa[::-1]
	
	def counts(self, seq, strand):
		aa = self.sequence(seq, strand)
		return Counter(aa)

	def frequencies(self, seq, strand):
		counts = self.counts(seq, strand)
		total = sum(counts.values())
		for aa in counts:
			counts[aa] = counts[aa] / total
		return counts
