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
    with open("/Users/wachsmannj/Documents/test_SKA2/cluster_15.txt") as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    if args["--fasta"]:
        print("run ska fasta")
        ska_cpp.run_ska(lines, 31)
    elif args["--align"]:
        print("run ska align")
        # ska_cpp.ska_align(args["file-list"])
        ska_cpp.ska_align(lines, 31)
    elif args["--map"]:
        print("run ska map")
        # ska_cpp.ska_map(args["file-list"])
    else:
        print("Option error!")
        sys.exit(1)

    # # list with all the contigs of all the isolates
    # seq_list = []
    # # list telling how many contigs where in each isolates
    # length_list = []
    #
    # # read in fasta files
    # input_path = '/Users/wachsmannj/Documents/test_SKA2/'
    # for file in os.listdir(input_path):
    #     global_file = '/'.join([input_path, file])
    #     count_records = 0
    #     # if re.search('\.fa$|\.fasta', str(file).strip()):
    #     if re.search('\.fasta$', str(file).strip()):
    #         # with open(global_file) as fasta_sequence:
    #         for record in SeqIO.parse(global_file, "fasta"):
    #             # print(record.id, record.seq)
    #             seq_list.append(str(record.seq))
    #             count_records += 1
    #         length_list.append(count_records)
    # ska_cpp.run_ska(seq_list, length_list, 7)
    # ska_cpp.run_ska(["ACTGAATC", "ACTCAATC", "ACTGAATC"], 7)

    sys.exit(0)


if __name__ == "__main__":
    main()
