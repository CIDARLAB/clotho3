from act_parser import act_parser
from act_query import act_query
from ClothoPy.protein_holder import Polypeptide
from uniprot_to_ncbi import uniprot_to_ncbi
from ClothoPy.accn_retrieval import call_accn
from TBlastN import TBlastN
from ClothoPy.blast_holder import Blast_Record
from protein_to_orf import _protein_to_orf
from design_operon import design_operon
from grab_rbs import grab_rbs
from ClothoPy.blast_retrieval import call_blast


"""
Converts a Pathway object into a list of ORFs.
"""
def pathway_to_orf(pathway):
	enzyme_list = pathway['reactions'][0][2]['enzymes']
	orf_list = []
	# print pathway['reactions']
	# for reaction in pathway['reactions']:
		# enzyme_list.append(chooseEnzyme(reaction))
	# what is being processed
	# for enzyme in enzyme_list:
	# 	for poly in enzyme['polypeptides']:
	#		print poly['uniprot']
	for enzyme in enzyme_list:
		#print enzyme
		orf = process_polypeptide(enzyme)
		if orf != [None]:
			orf_list += orf
	return orf_list

"""
For each Polypeptide in the Enzyme object, obtain its ORF.
"""
def process_polypeptide(enzyme):
	retriever = call_accn('nucleotide', 'gb', 'me@example.com')
	orfs = []
	for poly in enzyme['polypeptides']:
		# need to put in a check here.
		uni = uniprot_to_ncbi(poly['uniprot'])
		if uni != {}:
			#print "uniprot_to_ncbi"
			ncbi = uniprot_to_ncbi(str(poly['uniprot']))[str(poly['uniprot'])] # get ncbi number
			retriever.retrieve_gb([ncbi]) # get from ncbi
			record = retriever.records[0] # now we have the ncbi record
			orf = _protein_to_orf(record)
			orfs.append(orf)
			retriever.clearRecords()
		else:
			#print "tblastn"
			blaster = call_blast(poly['sequence'], 'me@example.com')
			record = blaster.retrieve_gb()
			# uni = BlastRecord(TBlastN([poly['sequence'], None]))
			# try:
			# 	ncbi = uni.alignments[0].accession # get ncbi number
			# 	retriever.retrieve_gb([ncbi]) # get from ncbi
			# 	record = retriever.records[0] # now we have the ncbi record
			# 	# need to add in data about start and end of the sequence in orf
			# 	# use query_start and query_end to get the sequence area
			if record is not None:
				orfs.append(record)
			#except Exception:
			#	print "No corresponding orf for %s" % poly.uniprot
	return orfs

"""
Trivial selector.
"""
def select_pathway(pathways):
	return pathways[0] #this is a trivial selection method, which should be changed

"""
Main function.
"""
def act_to_operon(query):
	#print "query"
	j = act_query(query)
	#print j
	#print "parser"
	paths = act_parser(j)
	#print paths
	#print "select"
	selected = select_pathway(paths)
	#print selected
	#print "orf list"
	orf_list = pathway_to_orf(selected)
	# for orf in orf_list:
	# 	print orf.GB.record
	#print "rbs list"
	rbs_list = grab_rbs()
	#print rbs_list
	#print "operon"
	operon = design_operon(orf_list, rbs_list, 6)
	return operon

def run(query):
	return act_to_operon(query)