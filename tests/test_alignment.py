import unittest
import ska_cpp
import numpy as np
from Bio import SeqIO

class MyTestCase(unittest.TestCase):
    def test_something(self):

        k = 31
        file1 = "/home/johanna/test_SKA2/integer_approach/MA_results/7553_5#73.skf"
        file2 = "/home/johanna/test_SKA2/integer_approach/MA_results/7622_2#17.skf"
        file3 = "/home/johanna/test_SKA2/integer_approach/MA_results/7553_5#55.skf"
        file4 = "/home/johanna/test_SKA2/integer_approach/MA_results/7553_5#93.skf"
        file5 = "/home/johanna/test_SKA2/integer_approach/MA_results/6925_2#81.skf"
        file6 = "/home/johanna/test_SKA2/integer_approach/MA_results/7622_5#64.skf"

        list1 = ska_cpp.compare_kmers(file1, k)
        list2 = ska_cpp.compare_kmers(file2, k)
        list3 = ska_cpp.compare_kmers(file3, k)
        list4 = ska_cpp.compare_kmers(file4, k)
        list5 = ska_cpp.compare_kmers(file5, k)
        list6 = ska_cpp.compare_kmers(file6, k)

        lists = [list1, list2, list3, list4, list5, list6]

        dict1 = {}
        dict2 = {}
        dict3 = {}
        dict4 = {}
        dict5 = {}
        dict6 = {}

        dicts = [dict1, dict2, dict3, dict4, dict5, dict6]

        for i in range(0, len(lists)):
            make_dicts(lists[i], dicts[i])

        merge_list = []

        for this_dict in dicts:
            for kmer in this_dict:
                merge_list.append(kmer)

        merge_list = set(merge_list)

        alignemnt_kmers = []
        for kmer in merge_list:
            count = 0
            for mydict in dicts:
                if kmer in mydict:
                    count += 1
            if (count/6 >= 0.9):
                alignemnt_kmers.append(kmer)

        alignment = []

        for i in range(0, len(alignemnt_kmers)): # kmer

            row = []
            for j in range(0, len(dicts)): #mydict
                if alignemnt_kmers[i] in dicts[j]:
                    row.append(dicts[j][alignemnt_kmers[i]])
                else:
                    row.append("-")
            if row.count(row[0]) != len(row):
                alignment.append(row)

        alignment1 = ([[row[i] for row in alignment] for i in range(len(alignment[0]))])

        alignment_file = []

        with open("/home/johanna/test_SKA2/integer_approach/MA_results/cluster28.aln") as handle:

            for record in SeqIO.parse(handle, "fasta"):
                print(record.seq)
                msa = []
                for base in record.seq:
                    msa.append(base)
                alignment_file.append(msa)

        for i in range(0, len(alignment1)):
            print("Test", i)
            self.assertEqual(sorted(alignment1[i]), sorted(alignment_file[i]))  # add assertion here


def make_dicts(list, mydict):
    for both in list:
        kmer = both.split(":")[0]
        base = both.split(":")[1]
        mydict[kmer] = base
    return mydict

if __name__ == '__main__':
    unittest.main()

