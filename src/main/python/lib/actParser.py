# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# For more information the Genbank file format, see:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245039/
# https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
# ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt

import json
from actQuery import actQuery

chem_dict = {}

"""
This is the main function, which takes in a dictionary of the result from act.20.
"""
def actParser(j):
	global chem_dict
	final_pathways = []
	pathways = j['results'][0]['pathways']['pathways']
	chemicals = j['results'][0]['pathways']['chemicals']

	description = j['results'][0]['pathways']['description']
	target = j['results'][0]['pathways']['target']
	version = j['results'][0]['pathways']['version']
	name = target + " Single Pathway"
	
	for chemical in chemicals:
		if 'SMILES' in chemical.keys():
			# in both substrates and products, they're strings, so I store them as strings as a result
			# change float to integer to string of integer
			chem_dict[str(int(chemical['id']))] = ('SMILES', chemical['SMILES'], chemical['name'])
		elif 'InChI' in chemical.keys():
			# change float to integer to string of integer
			chem_dict[str(int(chemical['id']))] = ('InChI', chemical['InChI'], chemical['name'])
	for path in pathways:
		single_pathway = {'description': description, 'target': target, 'version': version, 'name': name}
		interms = []
		reacts = []
		count = 0
		reactions = path['rxns_steps']
		for reaction in reactions:
			intermediates = reaction['substrates']
			products = reaction['products']
			reaction_desc = reaction['readable']

			# add this reaction's substrates to interms
			interm_temp, count = makeIntermediates(interms, intermediates, count)
			interms = interms + interm_temp

			# add this reaction's products to interms
			react_temp, count = makeIntermediates(interms, products, count)
			interms = interms + react_temp

			#make reactions to add to reacts
			r = makeReaction(interm_temp, react_temp, reaction['sequences'], reaction_desc)

			for i in range(len(interm_temp)):
				for j in range(len(react_temp)):
					reacts.append( [ interm_temp[i]['index'], react_temp[j]['index'], r ] )
		single_pathway['intermediates'] = interms
		single_pathway['reactions'] = reacts
		single_pathway['schema'] = "org.clothocad.model.SinglePathway"
		final_pathways.append(single_pathway)
	return final_pathways

"""
Create intermediates.
"""
def makeIntermediates(interms, intermediates, count):
	global chem_dict

	interm_temp = []

	# if intermediates['chemicals'] happens to be a single string
	if isinstance(intermediates['chemicals'], basestring):
		chem = intermediates['chemicals']
		interm = {"index": count}
		if chem in chem_dict.keys():
			if not inIntermediates(interms, chem_dict[chem][1]) and not inIntermediates(interm_temp, chem_dict[chem][1]):
				interm[chem_dict[chem][0]] = chem_dict[chem][1]		# InChI or SMILES : 
				interm['name'] = chem_dict[chem][2]					# name :
				interm_temp.append(interm)
				count += 1
		else:
			if not inIntermediates2(interms, chem) and not inIntermediates(interm_temp, chem):
				interm['name'] = chem
				interm_temp.append(interm)
				count += 1
	# usually intermediates['chemicals'] should be an array of chemicals
	else:
		for chem in intermediates['chemicals']:
			interm = {"index": count}
			if chem in chem_dict.keys():
				if not inIntermediates(interms, chem_dict[chem][1]) and not inIntermediates(interm_temp, chem_dict[chem][1]):
					interm[chem_dict[chem][0]] = chem_dict[chem][1]		# InChI or SMILES : 
					interm['name'] = chem_dict[chem][2]					# name :
					interm_temp.append(interm)
					count += 1
			else:
				if not inIntermediates2(interms, chem) and not inIntermediates(interm_temp, chem):
					interm['name'] = chem
					interm_temp.append(interm)
					count += 1

	# there isn't always a 'cofactors' field
	if 'cofactors' in intermediates.keys():
		# if intermediates['cofactors'] happens to be a single string
		if isinstance(intermediates['cofactors'], basestring):
			chem = intermediates['cofactors']
			interm = {"index": count}
			if chem in chem_dict.keys():
				if not inIntermediates(interms, chem_dict[chem][1]) and not inIntermediates(interm_temp, chem_dict[chem][1]):
					interm[chem_dict[chem][0]] = chem_dict[chem][1]		# InChI or SMILES : 
					interm['name'] = chem_dict[chem][2]					# name :
					interm_temp.append(interm)
					count += 1
			else:
				if not inIntermediates2(interms, chem) and not inIntermediates(interm_temp, chem):
					interm['name'] = chem
					interm_temp.append(interm)
					count += 1
		# usually intermediates['chemicals'] should be an array of chemicals
		else:
			for chem in intermediates['cofactors']:
				interm = {"index": count}
				if chem in chem_dict.keys():
					if not inIntermediates(interms, chem_dict[chem][1]) and not inIntermediates(interm_temp, chem_dict[chem][1]):
						interm[chem_dict[chem][0]] = chem_dict[chem][1]		# InChI or SMILES : 
						interm['name'] = chem_dict[chem][2]					# name :
						interm_temp.append(interm)
						count += 1
				else:
					if not inIntermediates2(interms, chem) and not inIntermediates(interm_temp, chem):
						interm['name'] = chem
						interm_temp.append(interm)
						count += 1
	return (interm_temp, count)

"""
Is this in intermediates?
Specifically for InChI and SMILES check.
"""
def inIntermediates(intermediates, check):
	for i in intermediates:
		if 'InChI' in i.keys():
			if i['InChI'] == check:
				return True
		if 'SMILES' in i.keys():
			if i['SMILES'] == check:
				return True
	return False

"""
Is this in intermediates?
Specifically for name check, aka no InChI or SMILES.
"""
def inIntermediates2(intermediates, name):
	for i in intermediates:
		if i['name'] == name:
			return True
	return False

"""
Make a reaction.
"""
def makeReaction(reactants, products, enzymes, description):
	fin = {"schema": "org.clothocad.model.Reaction", "description": description}
	reacts = [] #create the reactants
	for react in reactants:
		if 'InChI' in react.keys():
			reacts.append( "InChI=" + react['InChI'] )
		if 'SMILES' in react.keys():
			reacts.append( "SMILES=" + react['SMILES'] )
		else:
			reacts.append( react['name'] )
	fin['reactants'] = reacts
	prods = [] #create the products
	for prod in products:
		if 'InChI' in prod.keys():
			prods.append( "InChI=" + prod['InChI'] )
		if 'SMILES' in prod.keys():
			prods.append( "SMILES=" + prod['SMILES'] )
	fin['products'] = prods
	fin['enzymes'] = makeEnzyme(enzymes)
	return fin

"""
Make the enzyme list for reaction.
"""
def makeEnzyme(enzymes):
	fin = []
	# sometimes enzymes is just a single string
	if isinstance(enzymes, basestring):
		temp = {"schema": "org.clothocad.model.Enzyme"}
		temp['enzyme'] = [makePolypeptide(enzymes)]
		temp['observations'] = []
		fin.append(temp)
	# but usually it is a list of strings
	else:
		for enzyme in enzymes:
			temp = {"schema": "org.clothocad.model.Enzyme"}
			temp['enzyme'] = [makePolypeptide(enzyme)]
			temp['observations'] = []
			fin.append(temp)
	return fin

"""
Make a polypeptide.
"""
def makePolypeptide(sequence): # should probably have some more inputs, like uniprot and such
	return {'polypeptides': [sequence], 'uniprot': '<uniprot id>', \
		'source': '<source>', 'schema': 'org.clothocad.model.Polypeptide'}