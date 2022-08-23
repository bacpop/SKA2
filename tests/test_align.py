import unittest
import ska_cpp
from Bio import SeqIO

from ska.__main__ import read_in_files

class MyTestCase(unittest.TestCase):

    def test_ska_align(self):
        # paths, names = [], []
        # with open("10_test_cluster.txt", newline='\n') as file:
        #     lines = file.readlines()
        #     print(lines)
        #     for line in lines:
        #         paths.append(line.split("\t")[1].rstrip())
        #         names.append(line.split("\t")[0].rstrip())
        #         print(paths)
        #         print(names)
        paths, names = read_in_files("10_test_cluster.txt")
        ska_cpp.run_ska(paths, names, 3)

        ska_seqs = []
        with open("variant_alignment.aln") as fasta_sequences:
            for fasta in SeqIO.parse(fasta_sequences, "fasta"):
                sequence = str(fasta.seq)
                ska_seqs.append(sequence)

        actual_seqs = ["GCTTGTG", "GCTTT-G", "C-TTGTG"]

        for i in range(0, len(ska_seqs)):
            self.assertEqual(''.join(sorted(ska_seqs[i])), ''.join(sorted(actual_seqs[i])))  # add assertion here


if __name__ == '__main__':
    unittest.main()
