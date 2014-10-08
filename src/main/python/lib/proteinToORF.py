from Bio import Entrez, SeqIO
from ClothoPy.ProteinRetrieval import CallAccn
from ClothoPy.ProteinHolder import Polypeptide

def _ProteinToORF(prot):

    taken = prot.record
    coded_by = []
    trigger = False

    for feat in taken.features:
        if feat.type == 'CDS':
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
            records.append(Polypeptide(seq_record))
        result_handle.close()
    except:
        return None

    return records[0]

def run(*terms):
    return map(_ProteinToORF, terms)