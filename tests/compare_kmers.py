import unittest
import ska_cpp
import pandas as pd
import collections
import struct
import numpy as np


class MyTestCase(unittest.TestCase):

    def test_something(self):

        print("binary")
        f = open("/home/johanna/test_SKA2/integer_approach/new_test_cluster1/compare_k_mers/6925_2#80.skf", "r")
        a = np.fromfile(f, dtype=np.uint64)
        print("binary", a)
        f.close()

        k = 31

        file_ska = "/home/johanna/test_SKA2/integer_approach/new_test_cluster1/compare_k_mers/ska_pass_human.tsv"
        file_ska2 = "/home/johanna/test_SKA2/integer_approach/new_test_cluster1/compare_k_mers/6925_2#80.skf"

        ska2_kmer_list = sorted(ska_cpp.compare_kmers(file_ska2, k))[0:10]
        kmers_df = pd.read_csv(file_ska, sep='\t')
        ska_kmer_list = sorted(kmers_df[kmers_df.columns[0]].values.tolist())[0:10]
        # ska_kmer_list = sorted(ska_kmer_list)
        # ska2_kmer_list = sorted(ska2_kmer_list)
        kmerlist = []
        for kmer in ska2_kmer_list:
            kmer = kmer.split(":")[0]
            kmerlist.append(kmer)
        self.assertEqual(ska_kmer_list, kmerlist)  # add assertion here
        print("Test 1 successfully passed!")

        file_ska = "/home/johanna/test_SKA2/integer_approach/new_test_cluster1/compare_k_mers/ska_fail_human.tsv"
        file_ska2 = "/home/johanna/test_SKA2/integer_approach/new_test_cluster1/compare_k_mers/6925_1#53.skf"

        ska2_kmer_list = ska_cpp.compare_kmers(file_ska2, k)
        kmers_df = pd.read_csv(file_ska, sep='\t')
        ska_kmer_list = kmers_df[kmers_df.columns[0]].values.tolist()
        ska_base_list = kmers_df[kmers_df.columns[1]].values.tolist()
        ska_both = []

        ska_dict = {}
        for i in range(0, len(ska_base_list)):
            ska_dict[ska_kmer_list[i]] = ska_base_list[i]

        self.assertEqual(sorted(ska_kmer_list), sorted(ska2_kmer_list))
        print("Test 2 successfully passed!")

        ambiguity_bases = {
            "A": ["A"],
            "C": ["C"],
            "G": ["G"],
            "T": ["T"],
            "M": ["A", "C"],
            "R": ["A", "G"],
            "W": ["A", "T"],
            "S": ["C", "G"],
            "Y": ["T", "C"],
            "K": ["T", "G"],
            "V": ["A", "C", "G"],
            "H": ["A", "C", "T"],
            "D": ["A", "T", "G"],
            "B": ["T", "C", "G"],
            "N": ["A", "C", "G", "T"]
        }

        print("Ns", ska_base_list.count('N'))
        # ska_dict from ska
        for kmer in ska2_kmer_list:
            ska2_kmer = kmer.split(":")[0]
            ska2_base = kmer.split(":")[1]
            ska1_base = ska_dict[ska2_kmer]
            if ska1_base in ambiguity_bases[ska2_base] and ska1_base != 'N':
                match = True
            else:
                print(ska1_base)
                if ska1_base != 'N':
                    match = False
                    print(ska1_base, ska2_base, ambiguity_bases[ska2_base])
                    print("ska kmer", kmer)
                    print("ska_dict[ska2_kmer]", ska_dict[ska2_kmer])
            self.assertEqual(match, True)


if __name__ == '__main__':
    unittest.main()
