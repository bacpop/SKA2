from Bio.Seq import Seq
import unittest
import ska_cpp
import math
from Bio import SeqIO


class MyTestCase(unittest.TestCase):
    # def test_1_2(self):
    #     # Test 1: test binary conversion and reverse complement
    #     bitstring = ska_cpp.to_binary("ACCATTGCACA", 11)
    #     self.(343620, bitstring)
    #
    #     bitstring = ska_cpp.to_binary("ACCATTGCACACGTAAAACTT", 21)
    #     self.assertEqual(360312127519, bitstring)
    #
    #     bitstring = ska_cpp.to_binary("ACCATTGCACACGTAAAACTTTTTCCCGTAA", 31)
    #     self.assertEqual(377814649426400688, bitstring)
    #
    #     bitstring = ska_cpp.to_binary("ACTGCTA", 7)
    #     self.assertEqual(1948, bitstring)
    #
    #     # Test 2: test reverse complement
    #     return_smaller = ska_cpp.ReverseComp64(bitstring, 7)
    #     self.assertEqual(1948, return_smaller)
    #     # 00110010
    #     print("Test 1 & 2 successful!")
    #
    # def test_3_4_ambiguous_bases(self):
    #     # Test 3&4: test ambiguous_bases and lookuptable and lookuptable2
    #     split_bases = ['A', 'C', 'G', 'T']
    #     old_bases = ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']
    #     ambiguity_vector_calculated = []
    #     for split_base in split_bases:
    #         for old_base in old_bases:
    #             ambiguity_base = ska_cpp.test_ambiguity_bases(old_base, split_base)
    #             ambiguity_vector_calculated.append(ambiguity_base)
    #
    #     vector_ambiguity_truth = ['A', 'M', 'R', 'W', 'M', 'R', 'W', 'V', 'H', 'D', 'V', 'H', 'D', 'N', 'N',
    #                               'M', 'C', 'S', 'Y', 'M', 'V', 'H', 'S', 'Y', 'B', 'V', 'H', 'N', 'B', 'N',
    #                               'R', 'S', 'G', 'K', 'V', 'R', 'D', 'S', 'B', 'K', 'V', 'N', 'D', 'B', 'N',
    #                               'W', 'Y', 'K', 'T', 'H', 'D', 'W', 'B', 'Y', 'K', 'N', 'H', 'D', 'B', 'N']
    #     self.assertEqual(vector_ambiguity_truth, ambiguity_vector_calculated)
    #     print("Test 3 & 4 successful!")
    #
    # def test_5_check_N(self):
    #     # Test 5: test if Ns are correctly skipped
    #     k_test5 = 3
    #     skip_n_1 = ska_cpp.test_check_for_Ns("AAANAACT", k_test5)
    #     results_skip_n_1 = ["AAA", "AAC", "ACT"]
    #     self.assertEqual(results_skip_n_1, skip_n_1)
    #
    #     skip_n_2 = ska_cpp.test_check_for_Ns("AAANNNNAACT", k_test5)
    #     results_skip_n_2 = ["AAA", "AAC", "ACT"]
    #     self.assertEqual(results_skip_n_2, skip_n_2)
    #
    #     skip_n_3 = ska_cpp.test_check_for_Ns("AGGNNNACCCANNNNAACT", k_test5)
    #     results_skip_n_3 = ["AGG", "ACC", "CCC", "CCA", "AAC", "ACT"]
    #     self.assertEqual(results_skip_n_3, skip_n_3)
    #
    #     skip_n_4 = ska_cpp.test_check_for_Ns("AAGGGGGNACNCCNANACT", k_test5)
    #     results_skip_n_4 = ["AAG", "AGG", "GGG", "GGG", "GGG", "ACT"]
    #     self.assertEqual(results_skip_n_4, skip_n_4)
    #
    #     check = ska_cpp.check_for_N("ACTGCTA", 0)
    #     self.assertEqual([0, 0], check)
    #
    #     check = ska_cpp.check_for_N("AGGGNATAGCTAGCTAGAGATCTT", 0)
    #     self.assertEqual([1, 5], check)
    #     print("Test 5 successful!")
    #
    # #       TODO: check that in original code the position returned is supposed to be the one where the next k-mer starts
    #
    # def test_6_rolling_kmers(self):
    #     # Test 6: checking if the k-mers are created correctly
    #     k = 7
    #     sequence1 = "TATATNATTATATTATANNNNNNNN"
    #     c_dict = ska_cpp.test_new_rolling_kmer_binary(sequence1, k).keys()
    #     python_dict = rolling_kmers(sequence1, k).keys()
    #     self.assertEqual(c_dict, python_dict)
    #
    #     sequence2 = "TATATTCCCTTTTTCC"
    #     c_dict = ska_cpp.test_new_rolling_kmer_binary(sequence2, k).keys()
    #     python_dict = rolling_kmers(sequence2, k).keys()
    #     self.assertEqual(c_dict, python_dict)
    #
    #     sequence3 = "NNNNCANTATCTATATCTGAGGCGATCGAGGGGGCACCACACATATTAGAGGAGAGGAGAGAGAGTTTTCTCTCTCTCTCTGAGGGGGGGGGGGGACAANNNCACTANTA"
    #     c_dict = ska_cpp.test_new_rolling_kmer_binary(sequence3, k).keys()
    #     python_dict = rolling_kmers(sequence3, k).keys()
    #     self.assertEqual(c_dict, python_dict)
    #
    #     sequence4 = "GNNNACAANACATACATNNNCACTANTACNNATACATANNNA"
    #     c_dict = ska_cpp.test_new_rolling_kmer_binary(sequence4, k).keys()
    #     python_dict = rolling_kmers(sequence4, k).keys()
    #     self.assertEqual(c_dict, python_dict)
    #
    #     sequence5 = "NANAANACANAAAAANAATATATATATATATTATATAAAAAAAANAAA"
    #     c_dict = ska_cpp.test_new_rolling_kmer_binary(sequence5, k).keys()
    #     python_dict = rolling_kmers(sequence5, k).keys()
    #     self.assertEqual(c_dict, python_dict)
    #
    #     print("Test 6 successful!")

    def test_7_real_data(self):
        # Test 7: Tests if ska fasta works for small amounts of real data
        k = 31
        rfile = "/pp/PopPIPE_ska2/output/strains/28/rfile.txt"

        # read in from rfile
        paths, names = [], []
        with open(rfile, newline='\n') as file:
            lines = file.readlines()
            for line in lines:
                paths.append(line.split("\t")[1].rstrip())
                names.append(line.split("\t")[0].rstrip())

        # get input sequences into format to test: testing_run_ska(std::vector<std::string> seqs, int k)
        # seqs should be a list of sequences tha tall belong to the same isolate

        sequences = []
        for path in paths:
            current_sequence = []
            with open(path, 'r') as fasta_sequences:
                for fasta in SeqIO.parse(fasta_sequences, "fasta"):
                    current_contig = str(fasta.seq)
                    current_sequence.append(current_contig)
                sequences.append(current_sequence)

        # Run Python and C++ to return k-mers
        for isolate in sequences:
            python_dict = {}
            c_mers = []
            for seqs in isolate:
                python_dict = dict_rolling_kmers(seqs, k, python_dict)
            c_mers = ska_cpp.testing_run_ska(isolate, k, "kmer")

            python_mers = list(python_dict.keys())
            # print(sorted(c_mers))
            self.assertEqual(sorted(c_mers), sorted(python_mers))

        print("7.1 k-mer was successful!")

        # Run Python and C++ to return split base
        for isolate in sequences:
            python_dict = {}
            c_mers = []
            for seqs in isolate:
                python_dict = dict_rolling_kmers(seqs, k, python_dict)
            c_mers = ska_cpp.testing_run_ska(isolate, k, "both")

            items = []
            # ("item", list(python_dict.items()))

            for kmer in python_dict:
                item = kmer + ":" + python_dict[kmer]
                items.append(item)

            items = sorted(items)
            c_mers = sorted(c_mers)
            # print(python_dict)
            # print(c_mers)
            # python_mers = list(python_dict.values())
            # for something in range (0, len(items)):
            #     print(something, "python", items[something], "c++", c_mers[something], "\n")
            #     # self.assertEqual(items[something], c_mers[something])
            self.assertEqual(sorted(c_mers), sorted(items))

        print("Test 7.1 (both) successful!")

        # same thing as above just with a new input file:
        rfile = "/home/johanna/test_SKA2/integer_approach/MA_results/rfile.txt"
        rfile = "/pp/PopPIPE_ska2/output/strains/28/rfile.txt"

        # read in from rfile
        paths, names = [], []
        with open(rfile, newline='\n') as file:
            lines = file.readlines()
            for line in lines:
                paths.append(line.split("\t")[1].rstrip())
                names.append(line.split("\t")[0].rstrip())

        # get input sequences into format to test: testing_run_ska(std::vector<std::string> seqs, int k)
        # seqs should be a list of sequences tha tall belong to the same isolate

        sequences = []
        for path in paths:
            current_sequence = []
            with open(path, 'r') as fasta_sequences:
                for fasta in SeqIO.parse(fasta_sequences, "fasta"):
                    current_contig = str(fasta.seq)
                    current_sequence.append(current_contig)
                sequences.append(current_sequence)

        # Run Python and C++
        for isolate in sequences:
            python_dict = {}
            c_mers = []
            for seqs in isolate:
                python_dict = dict_rolling_kmers(seqs, k, python_dict)
            c_mers = ska_cpp.testing_run_ska(isolate, k, "kmer")

            python_mers = list(python_dict.keys())
            # upperlist = [x.upper() for x in sorted(python_mers)]

            self.assertEqual(sorted(c_mers), sorted(python_mers))

        print("Test 7.2 successful!")


look_up_table2 = ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N', 'X']

ambiguious_bases = [['A', 'M', 'R', 'W'],
                    ['M', 'C', 'S', 'Y'],
                    ['R', 'S', 'G', 'K'],
                    ['W', 'Y', 'K', 'T'],
                    ['M', 'M', 'V', 'H'],
                    ['R', 'V', 'R', 'D'],
                    ['W', 'H', 'D', 'W'],
                    ['V', 'S', 'S', 'B'],
                    ['H', 'Y', 'B', 'Y'],
                    ['D', 'B', 'K', 'K'],
                    ['V', 'V', 'V', 'N'],
                    ['H', 'H', 'N', 'H'],
                    ['D', 'N', 'D', 'D'],
                    ['N', 'B', 'B', 'B'],
                    ['N', 'N', 'N', 'N']]


def dict_rolling_kmers(sequence, k, current_kmer_dict):
    pos = 0
    current_kmer_dict = make_the_mers_ambiguity(pos, sequence, k, current_kmer_dict)
    return current_kmer_dict


def rolling_kmers(sequence, k):
    pos = 0
    current_kmer_dict = {}
    current_kmer_dict = make_the_mers(pos, sequence, k, current_kmer_dict)
    return current_kmer_dict


def reverse_compliment(inseq):
    # complement strand
    seq = inseq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    # reverse strand
    seq = seq[::-1]
    return seq


def init(pos, k, sequence):
    current_kmer = sequence[pos:pos + k]
    return current_kmer


def make_the_mers(pos, sequence, k, current_kmer_dict):
    # pos is the position of the beginning of the current k-mer
    # current k-mer
    # print("current sequence python:", sequence)
    current_kmer = init(pos, k, sequence)
    if 'N' in current_kmer:
        pos_of_N = current_kmer.index('N')
        pos = pos + pos_of_N + 1  # start with new k-mer directly after the N
        if pos < len(sequence) - k + 1:
            make_the_mers(pos, sequence, k, current_kmer_dict)
    else:
        split_kmer = current_kmer[0:math.floor(k / 2)] + current_kmer[math.floor(k / 2) + 1:k]
        split_base = current_kmer[math.floor(k / 2)]
        # add reverse compliment
        reverse_comp = reverse_compliment(split_kmer)
        if reverse_comp < split_kmer:
            # add split k-mer and split base to dictionary
            current_kmer_dict[reverse_comp] = reverse_compliment(split_base)
        else:
            # add split k-mer and split base to dictionary
            current_kmer_dict[split_kmer] = split_base
        # add ambiguity bases
        # TODO: ambiguity bases

        # add to position to get next k-mer
        pos += 1

        # getting next k-mer from new pos
        for i in range(pos, len(sequence) - k + 1):
            new_base = sequence[i + k - 1]
            if new_base == 'N':  # the new base is an N
                # update position
                pos = i + k
                if pos < len(sequence) - k + 1:
                    return make_the_mers(pos, sequence, k, current_kmer_dict)
                else:
                    return current_kmer_dict
            else:

                new_split_base = split_kmer[math.floor(k / 2)]
                s1 = split_kmer[1:math.floor(k / 2)]
                s2 = split_kmer[math.floor(k / 2) + 1:len(split_kmer)]
                split_kmer = s1 + split_base + s2 + new_base
                split_base = new_split_base

                # add reverse compliment
                reverse_comp = reverse_compliment(split_kmer)
                if reverse_comp < split_kmer:
                    reverse_base = reverse_compliment(split_base)
                    current_kmer_dict[reverse_comp] = reverse_base
                # add split k-mer and split base to dictionary
                else:
                    current_kmer_dict[split_kmer] = split_base
    return current_kmer_dict


def make_the_mers_ambiguity(pos, sequence, k, current_kmer_dict):
    # print(sequence)
    # pos is the position of the beginning of the current k-mer
    # current k-mer
    current_kmer = init(pos, k, sequence)
    if 'N' in current_kmer:
        pos_of_N = current_kmer.index('N')
        pos = pos + pos_of_N + 1  # start with new k-mer directly after the N
        if pos < len(sequence) - k + 1:
            make_the_mers_ambiguity(pos, sequence, k, current_kmer_dict)
    else:
        split_kmer = current_kmer[0:math.floor(k / 2)] + current_kmer[math.floor(k / 2) + 1:k]
        split_base = current_kmer[math.floor(k / 2)]
        # add reverse compliment

        reverse_comp = reverse_compliment(split_kmer)
        # reverse compliment
        if reverse_comp < split_kmer:
            # add split k-mer and split base to dictionary
            # ambiguity bases
            if reverse_comp in current_kmer_dict:
                x = look_up_table2.index((reverse_compliment(split_base)))
                y = look_up_table2.index(current_kmer_dict[reverse_comp])
                current_kmer_dict[reverse_comp] = ambiguious_bases[y][x]
                # print(x, "+", y, "=", ambiguious_bases[y][x])
                # print("1", split_kmer, ":", split_base, "saved as", reverse_comp, ":", ambiguious_bases[y][x])
            else:
                current_kmer_dict[reverse_comp] = reverse_compliment(split_base)
                # print("2", split_kmer, ":", split_base, "saved as", reverse_comp, ":", reverse_compliment(split_base))
        else:
            # add split k-mer and split base to dictionary
            if split_kmer in current_kmer_dict:
                x = look_up_table2.index(split_base)
                y = look_up_table2.index(current_kmer_dict[split_kmer])
                current_kmer_dict[split_kmer] = ambiguious_bases[y][x]
                # print(x, "+", y, "=", ambiguious_bases[y][x])
                # print("3", split_kmer, ":", split_base, "saved as", split_kmer, ":", ambiguious_bases[y][x])
            else:
                current_kmer_dict[split_kmer] = split_base
                # print("4", split_kmer, ":", split_base, "saved as", split_kmer, ":", split_base)

        # add to position to get next k-mer
        pos += 1

        # getting next k-mer from new pos
        for i in range(pos, len(sequence) - k + 1):
            new_base = sequence[i + k - 1]
            if new_base == 'N':  # the new base is an N
                # update position
                pos = i + k
                if pos < len(sequence) - k + 1:
                    return make_the_mers_ambiguity(pos, sequence, k, current_kmer_dict)
                else:
                    return current_kmer_dict
            else:

                new_split_base = split_kmer[math.floor(k / 2)]
                s1 = split_kmer[1:math.floor(k / 2)]
                s2 = split_kmer[math.floor(k / 2) + 1:len(split_kmer)]
                split_kmer = s1 + split_base + s2 + new_base
                split_base = new_split_base

                # add reverse compliment
                reverse_comp = reverse_compliment(split_kmer)

                if reverse_comp < split_kmer:
                    # add split k-mer and split base to dictionary
                    # ambiguity bases

                    if reverse_comp in current_kmer_dict:
                        # if split_kmer in current_kmer_dict:
                        x = look_up_table2.index((reverse_compliment(split_base)))
                        y = look_up_table2.index(current_kmer_dict[reverse_comp])
                        # print(x, "+", y, "=", ambiguious_bases[y][x])
                        current_kmer_dict[reverse_comp] = ambiguious_bases[y][x]
                        # print(x, "+", y, "=", ambiguious_bases[y][x])
                        # print("5", split_kmer, ":", split_base, "saved as", reverse_comp, ":", ambiguious_bases[y][x])
                    else:
                        current_kmer_dict[reverse_comp] = reverse_compliment(split_base)
                        # print("6", split_kmer, ":", split_base, "saved as", reverse_comp, ":", reverse_compliment(split_base))
                else:
                    # add split k-mer and split base to dictionary
                    if split_kmer in current_kmer_dict:
                        x = look_up_table2.index(split_base)
                        y = look_up_table2.index(current_kmer_dict[split_kmer])
                        current_kmer_dict[split_kmer] = ambiguious_bases[y][x]
                        # print(x, "+", y, "=", ambiguious_bases[y][x])
                        # print("7", split_kmer, ":", split_base, "saved as", split_kmer, ":", ambiguious_bases[y][x])
                    else:
                        current_kmer_dict[split_kmer] = split_base
                        # print("8", split_kmer, ":", split_base, "saved as", split_kmer, ":", split_base)

    return current_kmer_dict


if __name__ == '__main__':
    unittest.main()
