from Bio import Entrez, SeqIO
from ClothoPy.ProteinRetrieval import CallAccn
from ClothoPy.ProteinHolder import Polypeptide

def _ORFByName(terms):
    Entrez.email = 'nobody@example.com'
    call = CallAccn('protein', 'gb', 'nobody@example.com')
    term_string = terms[0] + '[ORGN] AND ' + terms[1] + '[PROT]'
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

    if len(call.records) == 0:
        return None

    taken = call.records[0].record
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
    key = spl[0]
    seqStart = spl[1].split(':')[1]
    seqEnd = spl[3]

    records = []

    try:
        result_handle = Entrez.efetch(db='nucleotide', rettype='gb', id=key, seq_start=seqStart, seq_stop=seqEnd)
        for seq_record in SeqIO.parse(result_handle, 'gb'):
            records.append(Polypeptide(seq_record))
            #print seq_record
        result_handle.close()
    except:
        return None

    return records[0]

def run(*terms):
    return map(_ORFByName, terms)