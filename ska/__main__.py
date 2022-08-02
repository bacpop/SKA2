"""ska: split k-mer alignment

Usage:
  __main__.py --fasta <fasta>
  __main__.py --align <align>
  __main__.py --map <map>

Options:
  -h --help         Show this help.
  --version         Show version.
  -kmer_size=<int>  Size of k-mers; needs to be off (default: 7)

"""
# -d --directory    Directory of fasta files (files have to end in *.fa or *.fasta)
#   ska fasta (--fasta) <files>... -o <output>
#   ska align -l <file-list> -o <output>
#   ska map ref -l <file-list> -o <output>
#   ska (-h | --help)
#   ska (--version)
#
# Options:
#   -h --help     Show this help.
#   --version     Show version.
#   ska_fasta     run ska fasta
#
#
#   -o <output>    Output prefix.
#   -l <file-list> File with a list of input files.
#
#
# """

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
#     print(arguments)
#     print("blub")
#     # TODO Check options here
#
#     return arguments

def main():
    args = docopt(__doc__, version = "ska version=" + __version__)
    # Create a database (sketch input)

    if args["--fasta"]:
        print("run ska fasta")
        # ska_cpp.run_ska(args["file-list"])
    elif args["--align"]:
        print("run ska align")
        # ska_cpp.ska_align(args["file-list"])
    elif args["--map"]:
        print("run ska map")
        # ska_cpp.ska_map(args["file-list"])
    else:
        print("Option error!")
        sys.exit(1)

    # list with all the contigs of all the isolates
    seq_list = []
    # list telling how many contigs where in each isolates
    length_list = []

    # read in fasta files
    input_path = '/Users/wachsmannj/Documents/test_SKA2/'
    for file in os.listdir(input_path):
        global_file = '/'.join([input_path, file])
        count_records = 0
        # if re.search('\.fa$|\.fasta', str(file).strip()):
        if re.search('\.fasta$', str(file).strip()):
            # with open(global_file) as fasta_sequence:
            for record in SeqIO.parse(global_file, "fasta"):
                # print(record.id, record.seq)
                seq_list.append(str(record.seq))
                count_records += 1
            length_list.append(count_records)
    ska_cpp.run_ska(seq_list, length_list, 7)
            # print(seq_list)

            # ska_cpp.run_ska(["ACTGAATC", "ACTCAATC", "ACTGAATC"], 7)

    sys.exit(0)


if __name__ == "__main__":
    main()
