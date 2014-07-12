# ID to Genbank

from ClothoPy.AccnRetrieval import CallAccn
import os

def _convertID(accession_id):
	retriever = CallAccn('nucleotide', 'gb', 'nobody@example.com')
	retriever.retrieve_gb([accession_id])
	#con = GBConverter(retriever.records[0])
	#con.convert()
	retriever.records[0].writeRecord('temp.gb')
	gen = open('temp.gb', 'rU').read()
	os.remove('temp.gb')
	return gen #con.d

def run(*accession_ids):
    return map(_convertID, accession_ids)
