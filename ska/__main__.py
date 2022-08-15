"""ska: split k-mer alignment

Usage:
  ska --fasta <fasta>
  ska --align <align>
  ska --map <map>

Options:
  -h --help         Show this help.
  --version         Show version.
  -kmer_size=<int>  Size of k-mers; needs to be off (default: 7)
  -o <output>       Output prefix.

"""


import os, sys, re
import ska_cpp
from Bio import SeqIO
from docopt import docopt

from .__init__ import __version__
#from __init__ import __version__


# def get_options():
#     from docopt import docopt
#     arguments = docopt(__doc__) #, version="ska version"+__version__)
#     # arguments = docopt(__doc__, version="ska fasta" + ska_cpp.run_ska)
#     # TODO Check options here
#
#     return arguments

def main():
    args = docopt(__doc__, version = "ska version=" + __version__)
    # Create a database (sketch input)
    paths, names = [], []
    with open("/Users/wachsmannj/Documents/test_SKA2/cluster_15.txt", newline='\n') as file:
        lines = file.readlines()
        for line in lines:
            paths.append(line.split("\t")[1].rstrip())
            names = [line.split("\t")[0].rstrip()]

    if args["--fasta"]:
        print("run ska fasta")
        ska_cpp.run_ska(paths, names, 31)
    elif args["--align"]:
        print("run ska align")
        # ska_cpp.ska_align(args["file-list"])
        # ska_cpp.ska_align(paths, names, 31)
    elif args["--map"]:
        print("run ska map")
        # ska_cpp.ska_map(args["file-list"])
    else:
        print("Option error!")
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
