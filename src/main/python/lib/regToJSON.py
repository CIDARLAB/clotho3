# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import xml.etree.ElementTree as ET

"""
Converts from registry (using the filename) to JSON.
"""
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

def regParse(inName, outName):
	root = regToJSON(inName)
	fin = {root.tag: [recTraverse(root)]}
	import json
	j = json.dumps(fin, indent=4)
	with open(outName, 'w') as f:
		f.write(j)

# regParse('registry/BBa_A340620.xml', 'BBa_A340620.json')
regParse('registry/BBa_B0000.xml', 'BBa_B0000.json')
regParse('registry/BBa_B0017.xml', 'BBa_B0017.json')
regParse('registry/BBa_B0104.xml', 'BBa_B0104.json')