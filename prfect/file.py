from genbank.file import File
from slippery.locus import Locus

class File(File):
	def construct_locus(self):
		return Locus()
