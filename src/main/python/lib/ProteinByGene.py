from Bio import Entrez, SeqIO
from ClothoPy.ProteinRetrieval import CallAccn

def _ProteinByGene(terms):
	Entrez.email = 'nobody@example.com'
    call = CallAccn('protein', 'gb', 'nobody@example.com')
    term_string = terms[0] + '[ORGN] AND ' + terms[1] + '[GENE]'
    ids = []
    lastlen = 0
    start = 0

    while True:
        handle = Entrez.esearch(db="protein", term=term_string, retstart=start, retmax=50)
        record = Entrez.read(handle)
        handle.close()
        ids = ids + record['IdList']
        if lastlen == len(ids):
            break
        else:
            lastlen = len(ids)
        start = start + 50

    call.retrieve_gb(ids)

	prot_list = "["
	for record in call.records:
		prot_list = "\n" + prot_list + record.toJSON() + ","
	prot_list = prot_list[:len(prot_list)-1] + "\n]"

	return prot_list

def run(*terms):
	return map(_ProteinByName, terms)