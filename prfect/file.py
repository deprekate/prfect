from genbank.file import File
from prfect.locus import Locus

class File(File, locus=Locus()):
	pass
