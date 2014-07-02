# Polynucleotide to Genbank

from ClothoPy.GenBankHolder import GenBank
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

def _convert_polynucleotide_to_genbank(poly):
	gen = GenBank()
	
	gen.description = poly.description
	gen.setName(poly.name)
	gen.setID(poly.id)
	gen.setDate(poly.submissionDate)
	gen.annotations['accessions'] = poly.accession
	gen.setVersion(int(poly.accession[len(poly.accession) - 1:len(poly.accession)]))
	
	for high in poly.highlights:
		gen.addFeatures(high.start, high.end, high.plusStrand, "misc_feature", "", high.description, \
			{'ApEinfo_fwdcolor':high.forColor, 'ApEinfo_revColor':high.revColor, \
			'inference':high.inference, 'label':high.description, 'note':high.notes}, \
			None, high.refSeq)
	
	typ = 'DNA'
	alpha = generic_dna
	if 'u' in poly.sequence or 'U' in poly.sequence:
		typ = 'RNA'
		alpha = generic_rna
	if 'm' in poly.sequence or 'M' in poly.sequence or 'r' in poly.sequence or 'R' in poly.sequence:
		typ = ''
		alpha = generic_protein
	
	gen.sequence = MutableSeq(poly.sequence, alpha)
	
	gen.firstLine = "LOCUS       " + gen.name + "          " + str(len(poly.sequence)) \
		+ " bp " + ("ss-" if poly.isSingleStranded else "  ") + typ + "    " + \
		("linear" if poly.isLinear else "circular") + "      " + poly.submissionDate

	for feat in gen.features:
		for key in feat.qualifiers.keys():
			if feat.qualifiers[key] == '' or feat.qualifiers[key] == None or feat.qualifiers[key] == []:
				del feat.qualifiers[key]

	gen.record = SeqRecord(gen.sequence, gen.id, gen.name, \
        gen.description, None, \
        gen.features, gen.annotations)
	return gen

def run(*accession_ids):
	return map(_convert_polynucleotide_to_genbank, accession_ids)