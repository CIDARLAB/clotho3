import urllib
import xml.etree.ElementTree as ET
#import time
#filename = "all_parts.fasta"
#id_array = []

def regToJSON(filename):
	tree = ET.parse(filename)
	root = tree.getroot()
	return root

def recTraverse(root):
	if len(root.getchildren()) == 0:
		return root.text
	else:
		tree = {}
		for child in root:
			if len(child.getchildren()) == 0:
				tree[child.tag] = recTraverse(child)
			elif child.tag not in tree.keys():
				tree[child.tag] = [recTraverse(child)]
			else:
				tree[child.tag].append(recTraverse(child))
		return tree

def regParse(inName):
	root = regToJSON(inName)
	fin = {root.tag: [recTraverse(root)], "schema":"org.registry.model.Part"}
	import json
	j = json.dumps(fin, indent=4)
	return j

def _grabRegistry(theID):
	base = "http://parts.igem.org/xml/part."
	url = base + theID
	filename = theID + ".xml"
	urllib.urlretrieve(url, filename)
	jReg = regParse(filename)
	return jReg
	
def run(*ids):
	return map(_grabRegistry, ids)

"""
with open(filename, 'r') as f:
	for line in f:
		split = line.split()
		if line.startswith(">BBa_"): #and split[0][1:] > 'BBa_K777107':
			first = split[0][1:]
			print first
			id_array.append(first)

for i in id_array:
	url = base + i
	print url
	urllib.urlretrieve(url, "registry/" + i + ".xml")
	time.sleep(.5)
"""