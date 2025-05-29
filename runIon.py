#Leo Patrick Mulholland 29.05.25 
#Queen's University Belfast
#Program to run many Cowan-Config-Average-Ionization calculations.

import json 
import argparse
from inputclass import Input
from runManyCowan import runManyCowan

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--json', help='Path of JSON for generating run.')
args = parser.parse_args()


def main(input:Input):
    
    runManyCowan(input)
    
    return 0 

if (not args.json):
    input_default = Input()
    default = json.dumps(input_default.__dict__,indent=1)
    print(default) 

else: 

    oo = args.json    
    with open(oo, 'r') as j:
        contents = json.loads(j.read())
    
    input = Input(**contents)
    main(input)