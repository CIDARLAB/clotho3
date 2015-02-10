from Bio import Entrez, SeqIO
from ClothoPy.genbank_holder import Genbank
from ClothoPy.protein_retrieval import call_accn

"""
Expecting a Polypeptide object.
"""
def _protein_to_orf(prot):

    retriever = call_accn('protein', 'gb', 'me@example.com')
    retriever.retrieve_gb([prot.id])
    taken = retriever.records[0].record
    coded_by = []
    trigger = False

    for feat in taken.features:
        if feat.type == 'CDS' and 'coded_by' in feat.qualifiers.keys():
            coded_by = feat.qualifiers['coded_by']
            trigger = True
        else:
            pass

    if not trigger:
        return None
    
    spl = coded_by[0].split('.')
    if '(' in coded_by[0]:
        key = spl[0].split('(')[1]
        seqStart = spl[1].split(':')[1]
        seqEnd = spl[3][:len(spl[3])-1]
    else:
        key = spl[0]
        seqStart = spl[1].split(':')[1]
        seqEnd = spl[3]

    records = []

    try:
        result_handle = Entrez.efetch(db='nucleotide', rettype='gb', id=key, seq_start=seqStart, seq_stop=seqEnd)
        for seq_record in SeqIO.parse(result_handle, 'gb'):
            records.append(Genbank(seq_record))
        result_handle.close()
    except:
        return None

    return records[0]

def run(term):
    return _protein_to_orf(term) #map(_protein_to_orf, terms)