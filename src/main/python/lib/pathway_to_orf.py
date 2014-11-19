# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# A set of functions used to convert from a SinglePathway into a list
# of ORFs as GenBank objects. Convert to "Nucleotide" by initializing
# NewGenBank objects with the .record of the GenBank objects.

from ClothoPy.ProteinHolder import Polypeptide
from pathways.uniprot_to_ncbi import uniprot_to_ncbi
from ClothoPy.ProteinRetrieval import CallAccn
from ClothoPy.BlastP import BlastP
from ClothoPy.BlastHolder import BlastRecord
from proteinToORF import _ProteinToORF

"""
Converts a Pathway object into a list of ORFs.
"""
def pathway_to_orf(pathway):
	enzyme_list = []
	orf_list = []
	for reaction in pathway.reactions:
		enzyme_list.append(chooseEnzyme(reaction))
	for enzyme in enzyme_list:
		orf = processPolypeptide(enzyme)
		if orf != [None]:
			orf_list += orf
	return orf_list

"""
Trivial selector to choose Enzyme.
"""
def chooseEnzyme(reaction):
	# for now, this is a trivial selector
	return reaction.enzymes[0]

"""
For each Polypeptide in the Enzyme object, obtain its ORF.
"""
def processPolypeptide(enzyme):
	retriever = CallAccn('protein', 'gb', 'me@example.com')
	orfs = []
	for polypeptide in enzyme.polypeptides:
		for poly in polypeptide:
			# need to put in a check here.
			uni = uniprot_to_ncbi(poly.uniprot)
			if uni != {}:
				ncbi = uniprot_to_ncbi(poly.uniprot)[poly.uniprot] # get ncbi number
				retriever.retrieve_gb([ncbi]) # get from ncbi
				record = retriever.records[0] # now we have the ncbi record
				orf = _ProteinToORF(record)
				orfs.append(orf)
				retriever.clearRecords()
			else:
				uni = BlastRecord(TBlastN([poly.sequence, None]))
				try:
					ncbi = uni.alignments[0].accession # get ncbi number
					retriever.retrieve_gb([ncbi]) # get from ncbi
					record = retriever.records[0] # now we have the ncbi record
					orfs.append(record)
				except Exception:
					print "No corresponding orf for %s" % poly.uniprot

	return orfs

def run(*pathways):
    return map(pathway_to_orf, pathways)