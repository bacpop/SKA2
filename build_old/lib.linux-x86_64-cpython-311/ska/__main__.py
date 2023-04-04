"""
SKA2: Split k-mer Analysis (Version 2.0.0)

Usage:
  ska fasta -i <in> -k <kmer> [-o <output> | -t <threads> | -c <value>]
  ska align -i <in> -k <kmer> [-o <output> | -t <threads> | -c <value>]
  ska map -i <in> -k <kmer> [-o <output> | -t <threads> | -c <value>]
  ska (-h | --help)
  ska --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -k --kmer         Length of kmers
  -i --input        Input assembly files
  -o --output       Output directory to save skf/alignment [default: .]
  -t --threads      Amount of threads to use
  --compression     Chose a compression technique for the skf [default: False]
"""

import sys, os
import ska_cpp
import subprocess
from docopt import docopt
from .__init__ import __version__
# import os, psutilpython setup.py install


def read_in_files(file_path):
    paths, names = [], []
    with open(file_path, newline='\n') as file:
        lines = file.readlines()
        for line in lines:
            paths.append(line.split("\t")[1].rstrip())
            names.append(line.split("\t")[0].rstrip())
    return paths, names


def main():
    args = docopt(__doc__, version="ska version = " + __version__)
    # print(args)
    k = int(args["<kmer>"])
    output_directory = args["<output>"]
    in_files = args["<in>"]
    compression = args["-c"]
    if compression:
        compression = [True, -1]
    else:
        compression = [False, 0]
    if args["fasta"]:
        print("run ska fasta with k-mer length", k)
        cluster = in_files.split('/')
        cluster = cluster[len(cluster)-1]
        cluster = cluster.split('.')[0]
        paths, names = read_in_files(in_files)
        ska_cpp.run_ska_fasta(paths, names, k, output_directory)
        if compression[0]:
            for i in names:
                file = output_directory + i + ".skf"
                if compression[1] == -1:
                    subprocess.run(["lz4", file, "-1"])
                else:
                    subprocess.run(["lz4", file, "-9"])
        # create file to continue with align
        # output_file = cluster + ".skf"
        f = open(os.path.join(os.path.dirname(in_files), cluster + ".skf"), "w")
        print(os.path.join(os.path.dirname(in_files), cluster + ".skf"))
        for i in range(0, len(paths)):
            path = os.path.join(output_directory, names[i] + ".skf")
            line = names[i] + "\t" + path + "\n"
            f.write(line)
            # print(names[i], path)
        f.close()
    elif args["align"]:
        print("run ska align with k-mer length", k)
        paths, names = read_in_files(in_files)
        cluster = os.path.basename(in_files).replace(".skf", "")
        ska_cpp.run_ska_align(paths, names, k, output_directory, cluster)
        print("Done!")


    elif args["map"]:
        print("run ska map")
        # ska_cpp.run_ska_map(paths, names, k, output_directory)
    else:
        print("Option error!")
        sys.exit(1)

    # process = psutil.Process(os.getpid())
    # print("Memory used:", process.memory_info().rss / 1024 ** 2)  # in MBs
    sys.exit(0)
if __name__ == "__main__":
    main()
