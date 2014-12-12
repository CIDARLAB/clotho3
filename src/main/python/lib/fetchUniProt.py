# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import urllib
import xml.etree.ElementTree as ET

"""
This function fetches a UniProt record using its ID. Only UniProt and not UniRef, etc.
"""
def fetchUniProt(ID):
	en = 'http://www.uniprot.org/uniprot/' #base url
	url = en + ID + ".xml"
	out = urllib.urlopen(url)
	s = out.read() #xml file
	out.close()

	root = ET.fromstring(s)
	entry = root.getchildren()[0]

	# these things you just have to go through uniprot xml files to realize
	sequence = "".join(entry.find('{http://uniprot.org/uniprot}sequence').text.split('\n'))
	organism = entry.find('{http://uniprot.org/uniprot}organism').getchildren()[0].text #scientific name

	protein = entry.find('{http://uniprot.org/uniprot}protein')
	recommendedName = protein.find('{http://uniprot.org/uniprot}recommendedName')
	name = recommendedName.find('{http://uniprot.org/uniprot}fullName').text

	return {'uniprot': ID, 'name': name, 'source': organism, 'sequence': sequence}