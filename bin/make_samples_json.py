#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names.
'''

import json
from glob import glob

# Change this line to match your filenames.
SRRS = [line.rstrip('\n') for line in open('meta_data/SRR_Acc_List.txt')]
FILES = {}

for srr_id in SRRS:
    FILES[srr_id] = {}
    FILES[srr_id]['R1'] = srr_id

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
