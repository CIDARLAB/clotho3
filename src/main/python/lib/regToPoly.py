import json
import random
from registryCollector import _grabRegistry

"""Assumes that registry is in JSON string."""
def _regToPoly(registry):
    fin = {}
    fin['isLinear'] = True
    fin['isSingleStranded'] = False
    if hasattr(registry, "part_short_desc"):
        fin['description'] = registry.part_short_desc
    if hasattr(registry, "part_name"):
        fin['accession'] = registry.part_name
        fin['name'] = registry.part_name
    if hasattr(registry, "part_entered"):
        MONTHS = {'01':'JAN', '02':'FEB', '03':'MAR', '04':'APR', '05':'MAY', '06':'JUN', \
        '07':'JUL', '08':'AUG', '09':'SEP', '10':'OCT', '11':'NOV', '12':'DEC'}
        dateSplit = change[key].split("-")
        date = dateSplit[2] + '-' + MONTHS[dateSplit[1]] + '-' + dateSplit[0]
        fin['submissionDate'] = date
    if hasattr(registry, "sequences"):
        fin['sequence'] = "".join(registry.sequences[0].seq_data.split())
        #I don't know how this will look like
    if hasattr(registry, "id"):
        fin['id'] = registry.id
    if hasattr(registry, "features"):
        fin['highlights'] = []
        #also don't know how this will be processed
        for f in registry.features[0].feature:
            r = lambda: random.randint(0,255)
            color1 = '#%02X%02X%02X' % (r(),r(),r())
            color2 = '#%02X%02X%02X' % (r(),r(),r())

            temp = {}
            temp['forColor'] = color1
            temp['recColor'] = color2
            temp['inference'] = None
            temp['description'] = None
            temp['notes'] = []
            temp['start'] = f['startpos']
            temp['end'] = f['endpos']
            if f['title'] == None or f['title'] == "":
                temp['name'] = "unnamed"
            else:
                temp['name'] = f['title']
            if f['direction'] == 'forward':
                temp['plusStrand'] = True
            else:
                temp['plusStrand'] = False
            if f['type'] != "BioBrick" or \
                (f['type'] == "BioBrick" and not f['title'].startswith('BBa_') \
                or f['title'] == change['part_name']):

                temp['refSeq'] = None
                fin['highlights'].append(temp)
            else:
                
                jj = _grabRegistry(f['title'])
                res = _regToPoly(jj)
                load = json.loads(res)
                concat = dict({'refSeq': load}.items() + temp.items())
                fin['highlights'].append(concat)
    return json.dumps(fin, indent=4)

def run(*ids):
    return map(_regToPoly, ids)