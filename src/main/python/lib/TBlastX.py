from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

def TBlastX(arr):
	result_handle = NCBIWWW.qblast("tblastx", "nt", Seq(arr[0]), alignments=arr[1], hitlist_size=arr[1])
	blast_record = NCBIXML.read(result_handle)
	return blast_record

def run(*arrs):
    return map(TBlastX, arrs)