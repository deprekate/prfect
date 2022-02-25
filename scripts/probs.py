import os
import sys
import random

from genbank.file import File

def is_hexa(seq):
	if 'n' in seq or 'N' in seq:
		return False
	elif (seq[0] == seq[1] == seq[2]) and (seq[3] == seq[4] == seq[5]):
		print(seq)
		return True
	return False

def gc(dna):
	dna = dna.upper()
	a = dna.count('A')
	c = dna.count('C')
	g = dna.count('G')
	t = dna.count('T')
	return (g+c) / (a+c+g+t)

bases = list("ACGT")

true = tot = 0
for i in range(1000):
	tot += 1
	seq = ""
	for j in range(6):
		seq += random.choice(bases)
	if is_hexa(seq):
		true += 1

print(true/tot)


filepath = sys.argv[1]
for filename in os.listdir(filepath):
	filename = "GCF_000872145.1_ViralProj27891_genomic.gbff.gz"
	filename = "GCF_000908955.1_ViralProj212947_genomic.gbff.gz"
	true = tot = 0
	dna = File(os.path.join(filepath, filename)).dna()
	for i in range(len(dna)-5):
		tot += 1
		if is_hexa(dna[i:i+6]):
			true += 1
	print(filename, gc(dna), true/tot, flush=True)
	exit()
