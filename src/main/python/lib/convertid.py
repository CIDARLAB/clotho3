# ID to Genbank

from ClothoPy.AccnRetrieval import CallAccn
from Bio import SeqIO
from StringIO import StringIO


def _convertID(accession_id):
	retriever = CallAccn('nucleotide', 'gb', 'nobody@example.com')
	retriever.retrieve_gb([accession_id])
	#con = GBConverter(retriever.records[0])
	#con.convert()
	
	#retriever.records[0].writeRecord('temp.gb')
	#gen = open('temp.gb', 'rU').read()
	#os.remove('temp.gb')
	out_handle = StringIO()
	SeqIO.write(retriever.records[0].record, out_handle, "gb")
	gb_data = out_handle.getvalue()

	return gb_data #con.d

def run(*accession_ids):
    return map(_convertID, accession_ids)
