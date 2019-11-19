# -*- coding: utf-8 -*-
"""
Description:
    For a given set of compounds representing the growth medium (or seeds). Create 2 reactions
    for each compounds to maintain consistency of the network for flux analysis.
    For each compounds create:

        An exchange reaction: this reaction consumes the compound in the
        compartment 'C-BOUNDARY' and produces the compound in the compartment 'e' extracellular

        A transport reaction: this reaction consumes the compound in the compartment
        'e' extracellular' and produces the compound in the compartment 'c' cytosol
        ex: for seed 'cpd-a'

    1/ check if cpd-a in padmetSpec, if not, copy from padmetRef.

    2/ create exchange reaction: ExchangeSeed_cpd-a_b: 1 cpd-a (C-BOUNDARAY) <=> 1 cpd-a (e)

    3/ create transport reaction: TransportSeed_cpd-a_e: 1 cpd-a (e) => 1 cpd-a (c)

    4/ create a new file if output not None, or overwrite padmetSpec

::

    usage:
        padmet_medium.py --padmetSpec=FILE
        padmet_medium.py --padmetSpec=FILE -r [--output=FILE] [-v]
        padmet_medium.py --padmetSpec=FILE --seeds=FILE [--padmetRef=FILE] [--output=FILE] [-v]

    options:
        -h --help     Show help.
        --padmetSpec=FILE    path to the padmet file to update
        --padmetRef=FILE    path to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
        --seeds=FILE    the path to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
        --output=FILE    If not None, pathname to the padmet file updated
        -r    Use to remove all medium from padmet
        -v   print info

"""
import docopt
from padmet.classes import PadmetSpec
from padmet.classes import PadmetRef
from padmet.utils.management import padmet_medium

def main():
    args = docopt.docopt(__doc__)
    if args["--seeds"]:
        with open(args["--seeds"], 'r') as f:
            seeds = [line.split("\t")[0] for line in f.read().splitlines()]
    else:
        seeds = None
    padmet = PadmetSpec(args["--padmetSpec"])
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    output = args["--output"]
    verbose = args["-v"]
    remove = args["-r"]

    if output is None:
        output = args["--padmetSpec"]

    if not remove and not seeds:
        g_m = padmet.get_growth_medium()
        print("List of growth medium:")
        if g_m:
            print(list(g_m))
        else:
            print("[]")
    else:
        padmet_medium.manage_medium(padmet, seeds, padmetRef, verbose)
        padmet.generateFile(output)


if __name__ == "__main__":
    main()

    
    