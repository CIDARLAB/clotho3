from ClothoPy.ProteinRetrieval import CallAccn

def _fetchProtein(polypeptide):
	return polypeptide.toJSON()

def run(*accession_ids):
	retriever = CallAccn('protein', 'gb', 'nobody@example.com')
    retriever.retrieve_gb(accession_ids)
    return map(_fetchProtein, retriever.records)