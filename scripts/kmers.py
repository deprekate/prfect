import sys


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


	

sets = dict()
seqs = list()
with open(sys.argv[1]) as fp:
	for line in fp:
		seqs.append(line.rstrip())
sets['six'] = find_motif(is_same, seqs, 6)
sets['hexa'] = find_motif(is_hexa, seqs, 6)
sets['five'] = find_motif(is_same, seqs, 5)
sets['twofour'] = find_motif(is_twofour, seqs, 6)
sets['four'] = find_motif(is_same, seqs, 4)
for key,value in sets.items():
	print(key, len(value))
print('other', len(seqs))


for i,seq in enumerate(seqs):
	#print(">seq", i)
	print(seq)

