from ClothoPy.NewGenBankHolder import NewGenBank
from ClothoPy.NewGBToJSON import NewGBConverter
from ClothoPy.GenBankHolder import GenBank

def _convert_file_to_nucseq(file_name):
	genbank = GenBank(file_name)
	#con = GBConverter(genbank)
	#con.convert()
	#genbank.writeRecord('temp.gb')
	#gen = open('temp.gb', 'rU').read()
	#os.remove('temp.gb')
	con = NewGBConverter(NewGenBank(gen.record))
    con.convert()
	return con.d

def run(*accession_ids):
    return map(_convert_genbank_to_nucseq, retriever.records)
