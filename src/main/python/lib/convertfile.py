# File to Genbank

from ClothoPy.GenBankHolder import GenBank
import os

def convertFile(file_name):
	genbank = GenBank(file_name)
	#con = GBConverter(genbank)
	#con.convert()
	genbank.writeRecord('temp.gb')
	gen = open('temp.gb', 'rU').read()
	os.remove('temp.gb')
	return gen #con.d
	# this returns the literal text inside of the record