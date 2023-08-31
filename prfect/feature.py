from genbank.feature import Feature
import textwrap


class Feature(Feature):
	def more(self):
		return "mooore"

	def end(self):
		if self.strand > 0:
			return self.right()
		else:
			return self.left()

	def nested_inside(self, other):
		if other and other.left() < self.left() and other.right() > self.right():
			return True
		return False

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

	def wwrite(self, outfile):
		outfile.write('     ')
		outfile.write( self.type.ljust(16) )
		if not self.strand > 0:
			outfile.write('coomplement(')
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
		for key,value in self.tags.items():
			for line in textwrap.wrap( '/' + str(key) + '=' + str(value) , 58):
				outfile.write('                     ')
				outfile.write(line)
				outfile.write('\n')





