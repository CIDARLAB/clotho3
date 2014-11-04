import urllib
import json

def actQuery(chemical):
	en = urllib.quote("".join(str({"chemical":chemical}).strip()))
	u = urllib.urlopen(('http://act.20n.com:27080/api/sampleFAB/_find?criteria=' + en).replace("%27", "%22"))
	s = u.read()
	u.close()
	j = json.loads(s)
	return j