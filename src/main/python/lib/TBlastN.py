# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# Using TBlastN.
# Input: array ([string sequence, int no_of_alignments] or [string sequence, null])

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from ClothoPy.blast_holder import Blast_Record
from Bio.Seq import Seq

def TBlastN(arr):
	result_handle = NCBIWWW.qblast("tblastn", "nt", Seq(arr[0]), alignments=arr[1], hitlist_size=arr[1])
	blast_record = NCBIXML.read(result_handle)
	return Blast_Record(blast_record, arr[0])

def run(arr):
    return TBlastN(arr) #map(TBlastN, arrs)