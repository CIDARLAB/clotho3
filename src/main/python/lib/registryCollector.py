import urllib
import xml.etree.ElementTree as ET

def regRoot(data):
	root = ET.fromstring(data)
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

def regParse(inString):
	root = regRoot(inString)
	fin = recTraverse(root)['part_list'][0]['part'][0]
	fin["schema"] = "org.registry.model.Part"
	fin['id'] = 'org.registry.part.' + fin['part_short_name']
	import json
	j = json.dumps(fin, indent=4)
	return j

def _grabRegistry(theID):
	base = "http://parts.igem.org/xml/part."
	url = base + theID
	filename = theID + ".xml"
	socket = urllib.urlopen(url)
	hold = ""
	for line in socket:
		hold += line
	jReg = regParse(hold)
	return jReg
	
def run(*ids):
	return map(_grabRegistry, ids)
