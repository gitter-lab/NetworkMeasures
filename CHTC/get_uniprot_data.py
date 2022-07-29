#!/usr/bin/env python

from time import sleep
import requests
import json
import sys
import os

if __name__ == '__main__':

    # uniprot url for node attributes
    url_base = 'https://rest.uniprot.org/uniprotkb/search?query=xref:string-'
    pid_file = sys.argv[1]
    results = {}

    with open(pid_file) as f:
        for pid in f:

            sleep(1)
            url = url_base + pid + '&format=json'
            result = requests.get(url).text
            result = json.loads(result)

            if len(result['results']) == 0:
                result = {}
            else:
                result = result['results'][0]

            results[pid] = result

    # save to file in directory
    with open(os.path.join('results', pid_file + '.json'), 'w') as f:
        json.dump(results, f)