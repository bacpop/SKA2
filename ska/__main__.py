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
import subprocess
from Bio import SeqIO
from docopt import docopt

from .__init__ import __version__


# def get_options():
#     from docopt import docopt
#     arguments = docopt(__doc__) #, version="ska version"+__version__)
#     # arguments = docopt(__doc__, version="ska fasta" + ska_cpp.run_ska)
#     # TODO Check options here
#     return arguments


def read_in_files(file_path):
    paths, names = [], []
    with open(file_path, newline='\n') as file:
        lines = file.readlines()
        for line in lines:
            paths.append(line.split("\t")[1].rstrip())
            names.append(line.split("\t")[0].rstrip())
    return paths, names


def main():
    args = docopt(__doc__, version="ska version=" + __version__)
    compression = [False, 0]
    # compression = [True, -1]
    #


    if args["--fasta"]:
        print("run ska fasta")
        paths, names = read_in_files("/Users/wachsmannj/Documents/test_SKA2/integer_approach/small_test/small_test_cluster.txt")
        # for i in range (0, 10):
        ska_cpp.run_ska_fasta(paths, names, 31)
        if compression[0]:
            for i in names:
                file = "/Users/wachsmannj/Documents/test_SKA2/integer_approach/testing/" + i + ".skf"
                if compression[1] == -1:
                    subprocess.run(["lz4", file, "-1"])
                else:
                    subprocess.run(["lz4", file, "-9"])
    elif args["--align"]:
        print("run ska align")
        paths, names = read_in_files("/Users/wachsmannj/Documents/test_SKA2/integer_approach/small_test/small_test_cluster_skf.txt")
        ska_cpp.run_ska_align(paths, names, 31)
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
