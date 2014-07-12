# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# For more information the Genbank file format, see:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245039/
# https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
# ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
#
# This prints out the NewGenBank object in a JSON format.

from NewGenBankHolder import NewGenBank
import json

class NewGBConverter:
    def __init__(self, gb):
        self.gb = gb #this is a NewGenBank object, not a GenBank object
        self.highs = []
        #self.annotes = {}
        self.d = {}
        self.json = ""

    def convert(self):
        self.highlight()
        self.d = {'schema': 'org.clothocad.model.Polynucleotide',
        'description': self.gb.description, \
        #'type': self.gb.type, \
        'sequence': self.gb.sequence, \
        'name': self.gb.name, \
        'id': self.gb.id, \
        'accession': self.gb.accn[0], \
        #'organism': self.gb.organism, \
        'date': self.gb.date, \
        'highlight': self.highs, \
        #'data_file_division': self.gb.data_file_division, \
        #'pubmed': self.gb.pubmed, \
        'isLinear': self.gb.isLinear, \
        'isSingleStranded': self.gb.isSingleStranded }

    def printJSON(self, indent):
        self.json = json.dumps(self.d, indent=indent)
        return self.json

    """
    def annotate(self):
        self.annotes = self.gb.annotations
        ref = self.gb.annotations['references']
        del self.annotes['references']
        self.annotes['references'] = []
        for r in ref:
            self.annotes['references'].append( \
            {'title': r.title, 'authors': r.authors, 'comment': r.comment, \
            'consrtm': r.consrtm, 'journal': r.journal, \
            'start': r.location[0].start.position, 'end': r.location[0].end.position, \
            'strand':r.location[0].strand, \
            'medline_id': r.medline_id, 'pubmed_id': r.pubmed_id} )
    """

    highlighted = {'ApEinfo_fwdcolor':'forColor', 'ApEinfo_revcolor':'revColor', \
    'inference':'inference', 'label': 'description', 'note': 'notes'}

    def highlight(self):
        high = self.gb.highlight
        count = 1;
        for f in high:
            loc = f.location
            qual = f.qualifiers
            self.highs.append( 
            {'start':loc.start.position, 'end': loc.end.position, \
            'forColor': '#000000', 'revColor': '#000000', 'inference': None, \
            'notes': None, 'description': None }) #notes is supposed to be the qualifier notes if it exists
            #'strand': loc.strand, 'type':f.type, \
            #'id:': f.id, \
            #'location operator':f.location_operator, 'ref': f.ref,\
            #'ref db': f.ref_db} )
            for q in qual.keys():
                if q in self.highlighted.keys():
                    self.highs[len(self.highs) - 1][self.highlighted[q]] = qual[q][0]
            count += 1

    """
    def letters(self):
        let = self.gb.letter_annotations
        for l in let:
            self.lets
    """
