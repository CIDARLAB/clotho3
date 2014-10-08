# Genbank to Polynucleotide (New Genbank)

from ClothoPy.GenBankHolder import GenBank
from ClothoPy.NewGenBankHolder import NewGenBank
from ClothoPy.NewGBToJSON import NewGBConverter
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO
import ClothoPy.ClothoSeqIO

def _convertGB(gb):
    #open('temp.gb', 'w').write(gb)
    #genbank = GenBank('temp.gb')
    #os.remove('temp.gb')

    gb_handle = StringIO(gb)
    record = ClothoPy.ClothoSeqIO.read(gb_handle, 'gb')
    
    #gb.record = SeqRecord(gb.sequence, id=gb.id, name=gb.name, \
    #   description=gb.description, dbxrefs=gb.dbxrefs, features=gb.features, \
    #   annotations=gb.annotations, letter_annotations=gb.letter_annotations)
    con = NewGBConverter(NewGenBank(record))
    con.convert()
    return con.d

def run(*accession_ids):
    return map(_convertGB, accession_ids)