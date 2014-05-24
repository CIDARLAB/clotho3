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
        self.feats = []
        self.annotes = {}
        self.d = {}
        self.json = ""

    def convert(self):
        self.d = {'description': self.gb.description, \
        'type': self.gb.type, \
        'sequence': self.gb.sequence, \
        'name': self.gb.name, \
        'ID': self.gb.id, \
        'organism': self.gb.organism, \
        'date': self.gb.date, \
        'data_file_division': self.gb.data_file_division, \
        'pubmed': self.gb.pubmed, \
        'isLinear': self.gb.isLinear, \
        'isSingleStranded': self.gb.isSingleStranded }

    def printJSON(self, indent):
        self.json = json.dumps(self.d, indent=indent)
        return self.json

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

    def feature(self):
        feat = self.gb.features
        count = 1;
        for f in feat:
            loc = f.location
            qual = f.qualifiers
            self.feats.append(
            {'start':loc.start.position, 'end': loc.end.position,\
            'strand': loc.strand, 'type':f.type, \
            'id:': f.id, \
            'location operator':f.location_operator, 'ref': f.ref,\
            'ref db': f.ref_db} )
            for q in qual.keys():
                self.feats[len(self.feats) - 1][q] = qual[q][0]
            count += 1

    def letters(self):
        let = self.gb.letter_annotations
        for l in let:
            self.lets
