# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

from find_local_hairpins import overall_scoring
from ClothoPy.genbank_holder import Genbank
from StringIO import StringIO
import ClothoPy.ClothoSeqIO
from Bio import Alphabet
from Bio.SeqRecord import SeqRecord

def poly_to_gen(poly):
    gen = Genbank()
    
    gen.description = ""
    gen.setName(poly["name"])
    gen.setID(poly["id"])
    date = ""
    gen.setDate(date)
    for high in poly["highlights"]:
        h_start = None
        h_end = None
        h_strand = None
        h_type = "misc_feature"
        h_loc = ""
        h_id = "<unknown id>"
        h_for = ""
        h_rev = ""
        h_inf = ""
        h_desc = ""
        h_note = []
        h_ref = ""

        h_start = high["start"]
        h_end = high["end"]
        h_strand = high["plusStrand"]
        h_for = high["forColor"]
        h_rev = high["revColor"]
        h_desc = high["description"]
        h_ref = high["refSeq"]

        gen.addFeatures(h_start, h_end, h_strand, "misc_feature", "", h_id, \
            {'ApEinfo_fwdcolor':h_for, 'ApEinfo_revColor':h_rev, \
            'inference':h_inf, 'label':h_desc, 'note':h_note}, \
            None, h_ref)
    
    gen.sequence = poly["sequence"]
    
    typ = ""
    if isinstance(gen.sequence.alphabet, Alphabet.DNAAlphabet):
        typ = "DNA"
    elif isinstance(gen.sequence.alphabet, Alphabet.RNAAlphabet):
        typ = "RNA"

    h_single = poly["isSingleStranded"]
    h_linear = poly["isLinear"]
    gen.firstLine = "LOCUS       " + gen.name + "          " + str(len(poly["sequence"])) \
        + " bp " + ("ss-" if h_single else "  ") + typ + "    " + \
        ("linear" if h_linear else "circular") + "      " + date

    for feat in gen.features:
        for key in feat.qualifiers.keys():
            if feat.qualifiers[key] == '' or feat.qualifiers[key] == None or feat.qualifiers[key] == []:
                del feat.qualifiers[key]

    gen.record = SeqRecord(gen.sequence, gen.id, gen.name, \
        gen.description, None, \
        gen.features, gen.annotations)

    out_handle = StringIO()
    ClothoPy.ClothoSeqIO.write(gen.record, out_handle, "gb")
    gb_data = out_handle.getvalue()

    return gb_data #con.d

def design_operon(orf_list, rbs_list, limit):
    operon = ""
    operon_name = ""
    highlights = []
    location = 0
    first = True

    for orf in orf_list:

        #print orf + " :: before" + str(len(orf))
        orf_seq = orf.sequence.lower()
        if len(orf_seq) >= 38:
            orf_seq = orf_seq[:37]
        #print orf + " :: after" + str(len(orf))
        min_score = float("infinity")
        min_rbs = ""

        for r in rbs_list:
            #print first
            is_intergenic_count = 0
            is_not_intergenic_count = 0
            if first:
                #print r
                if not r['isIntergenic']:
                    rbs = r['sequence']
                    temp_score = overall_scoring(rbs, orf_seq, limit, False)
                    if temp_score < min_score:
                        min_score = temp_score
                        min_rbs = r
                    is_not_intergenic_count += 1
            else:
                if r['isIntergenic']:
                    rbs = r['sequence']
                    temp_score = overall_scoring(rbs, orf_seq, limit, False)
                    if temp_score < min_score:
                        min_score = temp_score
                        min_rbs = r
                    is_intergenic_count += 1
        if first:
            first = False
        if is_intergenic_count > 1:
            rbs_list.remove(min_rbs) #remove entry from list of RBSs
        min_rbs_seq = min_rbs['sequence'].lower()
        if min_rbs_seq.endswith("atg") and orf_seq.startswith("atg"):
            min_rbs_seq = min_rbs_seq[:len(min_rbs_seq) - 3]
        addition = min_rbs_seq + orf_seq
        operon += addition
        #print min_rbs
        operon_name += min_rbs['name'] + "." + orf.name + "."

        import random
        r = lambda: random.randint(0,255)
        color1 = '#%02X%02X%02X' % (r(),r(),r())
        color2 = '#%02X%02X%02X' % (r(),r(),r())

        highlights.append( {'start': location, \
            'end': location + len(addition) - 1, \
            'refSeq': min_rbs['id'], \
            #'sequence': addition, \
            'schema': 'org.clothocad.model.Highlight', \
            'plusStrand': True, \
            'description': min_rbs['name'], \
            'forColor': color1, \
            'revColor': color2} )

        location += len(addition)

        # name should be a concatenation of rbs1.cds1.rbs2.cds2.etc
        # (orf_list is going to be a list of polynucleotides and not just a straight list)
        # isSingleStranded = False
        # isLinear = True
        # no accession number or submissionDate
        # description?

    return poly_to_gen ({ \
        'id': ".", \
        'sequence': operon, \
        'schema': 'org.clothocad.model.Polynucleotide', \
        'highlights': highlights, \
        'isLinear': True, \
        'isSingleStranded': False, \
        'name': operon_name[:len(operon_name)-1] })