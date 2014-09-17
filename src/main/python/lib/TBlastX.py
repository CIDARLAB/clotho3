from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from ClothoPy.BlastHolder import BlastRecord

def TBlastX(arr):
	result_handle = NCBIWWW.qblast("tblastx", "nt", Seq(arr[0]), alignments=arr[1], hitlist_size=arr[1])
	blast_record = NCBIXML.read(result_handle)
	return BlastRecord(blast_record, arr[0])

def run(*arrs):
    return map(TBlastX, arrs)