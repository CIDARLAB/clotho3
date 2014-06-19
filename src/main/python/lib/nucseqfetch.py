from ClothoPy.AccnRetrieval import CallAccn
from ClothoPy.NewGenBankHolder import NewGenBank
from ClothoPy.NewGBToJSON import NewGBConverter

def _convert_genbank_to_nucseq(genbank_obj):
    con = NewGBConverter(NewGenBank(genbank_obj.record))
    con.convert()
    return con.d

def run(*accession_ids):
    retriever = CallAccn('nucleotide', 'gb', 'nobody@example.com')
    retriever.retrieve_gb(accession_ids)
    return map(_convert_genbank_to_nucseq, retriever.records)
