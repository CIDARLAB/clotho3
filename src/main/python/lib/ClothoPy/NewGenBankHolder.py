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
from GenBankHolder import GenBank

class NewGenBank:
	def __init__ (self, record):
		self.GB = GenBank(record)

		self.name = self.GB.name #important
		self.description = self.GB.description #important
		self.id = self.GB.id #important

		self.sequence = self.GB.sequence.tostring() #important
		self.alphabet = self.GB.sequence.alphabet
		self.type = 'DNA' #important
		if isinstance(self.alphabet, RNAAlphabet):
			self.type = 'RNA'
		elif isinstance(self.alphabet, ProteinAlphabet):
			self.type = 'Protein'

		self.features = self.GB.features
		if 'references' in self.GB.annotations.keys():
			self.references = self.GB.annotations['references']
		self.pubmed = [] #important

		self.annotations = self.GB.annotations
		data_file_division = ''#important
		if 'data_file_division' in self.GB.annotations.keys():
			self.data_file_division = self.GB.annotations['data_file_division']
		self.organism = 'Unknown'
		if 'organism' in self.GB.annotations.keys(): #important
			self.organism = self.GB.annotations['organism']
		self.date = ''
		if 'date' in self.GB.annotations.keys(): #important
			self.date = self.GB.annotations['date']

		self.firstLine = self.GB.firstLine
		self.isLinear = True #important
		if 'circular' in self.firstLine.split():
			self.isLinear = False

		self.isSingleStranded = 'ss-' in self.firstLine

	"""Fairly self explanatory. Just extracts what the Pubmed IDs are,
	if there are any."""
	def extractRefs():
		for ref in self.references:
			pm = ref.pubmed_id
			if pm != "":
				self.pubmed.append(pm)
