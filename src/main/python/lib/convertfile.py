# File to Genbank

from ClothoPy.GenBankHolder import GenBank

def _convertFile(file_name):
	genbank = GenBank(file_name)
	#con = GBConverter(genbank)
	#con.convert()
	
	#genbank.writeRecord('temp.gb')
	#gen = open('temp.gb', 'rU').read()
	#os.remove('temp.gb')

	out_handle = StringIO()
	SeqIO.write(genbank.record, out_handle, "gb")
	gb_data = out_handle.getvalue()
	return gb_data 

	# this returns the literal text inside of the record

def run(*accession_ids):
    return map(_convertFile, accession_ids)