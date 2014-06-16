# Genbank to Polynucleotide (New Genbank)

from ClothoPy.GenBankHolder import GenBank
from ClothoPy.NewGenBankHolder import NewGenBank
from ClothoPy.NewGBToJSON import NewGBConverter

def convertGB(gb):
	open('temp.gb', 'w').write(gb)
	genbank = GenBank('temp.gb')
	os.remove('temp.gb')
	con = NewGBConverter(NewGenBank(genbank.record))
	con.convert()
    return con.d