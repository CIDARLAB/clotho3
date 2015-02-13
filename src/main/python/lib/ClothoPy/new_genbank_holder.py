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
# This class holds a new version of the Genbank file format to be used
# in NucSeq.

from Bio.Alphabet import ProteinAlphabet, DNAAlphabet, RNAAlphabet
from genbank_holder import Genbank

class New_Genbank:
	def __init__ (self, record):
		self.GB = Genbank(record) # A lot of the work is done already in GenBank
		
		self.name = self.GB.name
		self.description = self.GB.description

		self.sequence = str(self.GB.sequence)
	
		self.date = None
		if 'date' in self.GB.annotations.keys():
			self.date = self.GB.annotations['date']
		if self.GB.id is None or self.GB.id == "":
			self.id == self.GB.name
		else:
			self.id = self.GB.id
		if 'gi' in self.GB.annotations.keys():
			self.id = self.GB.annotations['gi']
		self.accn = None
		if 'accessions' in self.GB.annotations.keys():
			self.accn = self.GB.annotations['accessions']

		self.highlight = self.GB.features
		
		self.firstLine = self.GB.firstLine
		self.isLinear = True
		if 'circular' in self.firstLine.split():
			self.isLinear = False

		self.isSingleStranded = 'ss-' in self.firstLine

	"""Fairly self explanatory. Just extracts what the Pubmed IDs are,
	if there are any."""
	def extract_refs():
		for ref in self.references:
			pm = ref.pubmed_id
			if pm != "":
				self.pubmed.append(pm)

	def get_record():
		return self.GB.recordrecord