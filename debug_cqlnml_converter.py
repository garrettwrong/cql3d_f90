#!/usr/bin/env python
'''

In the process of running TRANSP, it can dump a debug nml.
At the time of writing this dump is enabled for MPI runs, to assist transition.
Developers may make this optional, perhaps only for debugging, in the future.

This namelist comes _from within cql3d_f90_ and so uses f90 derived type syntax.

The standalone xcql3d expects the legacy f77 syntax NML by request of Bob@CompX.

You can use this tool to convet to the legacy style.

Note, methods to read, write, and print both styles are implimented in the code.
Future developers may change to the new style if they wish.

You convert like so, assuming your dump is called debug_new:

./debug_cqlnml_converter.py debug_new cqlinput

I would also recomend that you set VERBOSE = 1, to re-enable debug prints.

Good luck.

'''

from __future__ import print_function

import sys


NML_SECTIONS = ['SETUP0%',
                'SETUP%',
                'TRSETUP%',
                'SOUSETUP%',
                'EQSETUP%',
                'RFSETUP%',
                'FRSETUP%']

def convert(in_fn, out_fn=None):
    if not out_fn:
        out_fn = in_fn

    with open(in_fn, 'r') as f:
        lines = f.readlines()

    res = []
    missing_freya = True
    for line in lines:
        # if we dont have freya section, we will need to add it later
        if 'frsetup' in line.lower():
            missing_freya = False
        # remove the _NML from section name
        line = line.replace('_NML','')

        # remove the dertype prefix
        for sec in NML_SECTIONS:
            # though we could just partition,
            # this check is slightly more robust, since
            #  we could conceivably have % in a txt field
            if sec.lower() in line.lower():
                _,_,line = line.partition('%')
                break # only going to be one nml section...
        res.append(line)
    if missing_freya:
        #add an empty freya block
        res.append('\n\n')
        res.append('&FRSETUP\n')
        res.append('&END\n\n')



    with open(out_fn, 'w') as f:
        f.write(''.join(res))

if __name__ == "__main__":

    out_fn = None
    if len(sys.argv)==3:
        out_fn = sys.argv[2]
    elif len(sys.argv) != 2:
        print("{x} input_filename <output_filename>".format(x=sys.argv[0]))
        exit(1)

    convert(sys.argv[1], out_fn)
