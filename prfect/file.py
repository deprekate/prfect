from genbank.file import File
from prfect.locus import Locus

class File(File):
	def construct_locus(self):
		return Locus()
