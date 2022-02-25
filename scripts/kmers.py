import sys

'''
	probabilities:
	six     = 0.00097 =  0.25**5
	fivetwo = 0.00097 = (0.25**4)*(0.25)
	five    = 0.0039  =  0.25**4
	hexa    = 0.0039  = (0.25**2)**2
	fourtwo = 0.0039  = (0.25**3)*(0.25)
	four    = 0.0156  =  0.25**3
'''

def best_motif(seqs, k):
	motifs = dict()
	for seq in seqs:
		for i in range(len(seq)-k):
			motifs[seq[i:i+k]] = motifs.get(seq[i:i+k], 0) + 1
	return max(motifs, key=motifs.get)

def is_hexa(seq):
	if (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
		return True

def is_twofour(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5]):
		return True
	elif (seq[0] == seq[1] == seq[2] == seq[3]) and (seq[4] == seq[5]):
		return True
	return False

def is_twofive(seq):
	if (seq[0] == seq[1] ) and (seq[2] == seq[3] == seq[4] == seq[5] == seq[6]):
		return True
	elif (seq[0] == seq[1] == seq[2] == seq[3] == seq[4]) and (seq[5] == seq[6]):
		return True
	return False

def is_repeated(seq):
	codons = dict()
	for i in range(0,len(seq),3):
		codon = seq[i:i+3]
		codons[codon] = True
	if len(codons.items()) < 2:
		return True
	else:
		return False

def is_same(seq):
	return seq == len(seq) * seq[0]

def has(motif, seq, k):
	for i in range(len(seq)-k+1):
		if motif(seq[i:i+k]):
			return True
	return False

def find_motif(motif, seqs, k):
	new = list()
	for i in range(len(seqs) - 1, -1, -1):
		if has(motif, seqs[i], k):
			new.append(seqs[i])
			del seqs[i]
	return new


	
kind = dict()
sets = dict()
seqs = list()
with open(sys.argv[1]) as fp:
	for line in fp:
		cols = line.rstrip().split('\t')
		seq = cols[2][20:]
		seqs.append(seq)
		kind[seq] = cols[1]
#--- done with seqs ---#
for seq in seqs:
	print(seq, end='\t')
	for func,k in [(is_same,6), (is_twofive,7),(is_hexa,6),(is_same,5),(is_twofour,6),(is_same,4)]:
		print(has(func, seq, k), end='\t')
	print()
exit()


sets['six'] = find_motif(is_same, seqs, 6)
sets['fivetwo'] = find_motif(is_twofive, seqs, 7)
sets['hexa'] = find_motif(is_hexa, seqs, 6)
sets['five'] = find_motif(is_same, seqs, 5)
sets['fourtwo'] = find_motif(is_twofour, seqs, 6)
sets['four'] = find_motif(is_same, seqs, 4)
sets['repeat'] = find_motif(is_repeated, seqs, 9)

for key,value in sets.items():
	counts = dict()
	for seq in value:
		counts[kind[seq]] = counts.get(kind[seq], 0) + 1
	print(key, len(value), counts)
print('other', len(seqs))

exit()
for i,seq in enumerate(seqs):
	print(">seq", i, sep='')
	print(seq)

