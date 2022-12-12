import unittest
import ska_cpp
from Bio import SeqIO
from ska.__main__ import read_in_files


class MyTestCase(unittest.TestCase):

    # def test_ska_align(self):
    #     # paths, names = [], []
    #     # with open("10_test_cluster.txt", newline='\n') as file:
    #     #     lines = file.readlines()
    #     #     print(lines)
    #     #     for line in lines:
    #     #         paths.append(line.split("\t")[1].rstrip())
    #     #         names.append(line.split("\t")[0].rstrip())
    #     #         print(paths)
    #     #         print(names)
    #     paths, names = read_in_files("10_test_cluster.txt")
    #     ska_cpp.testing_run_ska(paths, names, 3)
    #
    #     ska_seqs = []
    #     with open("variant_alignment.aln") as fasta_sequences:
    #         for fasta in SeqIO.parse(fasta_sequences, "fasta"):
    #             sequence = str(fasta.seq)
    #             ska_seqs.append(sequence)
    #
    #     actual_seqs = ["GCTTGTG", "GCTTT-G", "C-TTGTG"]
    #
    #     for i in range(0, len(ska_seqs)):
    #         self.assertEqual(''.join(sorted(ska_seqs[i])), ''.join(sorted(actual_seqs[i])))  # add assertion here

    #   test the alignemnt function to see if the correct middle bases are outputted
    def test_all_without_reverse_complement(self):

        sequence1 = ['AACAA', 'CCCCCCC', 'GGGGATATATATAG']
        sequence2 = ['AACAA', 'CGCCC', 'GGGAG']
        sequence3 = ['AACAA', 'CCTCC', 'GGGATATATAGGGGG']
        # sequence1 = ['AACAAA']
        # sequence2 = ['AACAAA']
        # sequence3 = ['AACATA']
        sequences = []

        sequences.append(sequence1)
        sequences.append(sequence2)
        sequences.append(sequence3)
        k = 3

        paths = ['/Users/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/seq1.fa',
                 '/Users''/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/seq2.fa',
                 '/Users/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/seq3.fa']
        paths2 = ['/Users/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/seq1.skf',
                 '/Users''/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/seq2.skf',
                 '/Users/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/seq3.skf']

        names = ["seq1", "seq2", "seq3"]

        # Python Code:
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

        # SKA FASTA
        vector_of_dicts = []
        for i in range(0, len(sequences)):
            print(sequences[i])
            current_dict = {}
            for j in range(0, len(sequences[i])):
                current_sequence = sequences[i][j]
                for position in range(0, len(current_sequence)-k+1):
                    # make k-mers
                    split_kmer = current_sequence[int(position): int(position+k/2)] + current_sequence[int(position+k/2+1): int(position+k)]
                    base = current_sequence[int(k/2 + position)]
                    print(split_kmer, base)

                    # check if kmers are already in current_dictionary
                    if split_kmer in current_dict:
                        x = look_up_table2.index(base)
                        y = look_up_table2.index(current_dict[split_kmer])
                        current_dict[split_kmer] = ambiguious_bases[y][x]
                    else:
                        current_dict[split_kmer] = base

                # add dict to vector to transfer
            vector_of_dicts.append(current_dict)

        print('vector of dicts', vector_of_dicts)
        # create unique_merge_vector to have all unique kmers to iterate through
        merge_vector = []
        for dic in vector_of_dicts:
            merge_vector.append(list(dic.keys()))
        flat_merge_vector = [item for blub in merge_vector for item in blub]
        myset = set(flat_merge_vector)
        unique_merge_vector = list(myset)
        unique_merge_vector.sort()

        # create alignment dictionary
        alignment_dict = {}

        # TODO: filter out:
        #     1. abundance above 0.5%
        #     2. if vector all the same bases (bc then there is no variance)
        for kmer in unique_merge_vector:
            abundance = 0
            variation = False
            previous = 'J'
            current_base_vector = []
            for isolate_dict in vector_of_dicts:
                # print('isolate_dict', isolate_dict)
                if kmer in isolate_dict:
                    # print('kmer', kmer)
                    abundance += 1
                    current_base_vector.append(isolate_dict[kmer])
                    if (previous != 'J' and previous != isolate_dict[kmer]):
                        variation = True
                    previous = isolate_dict[kmer]
                else:
                    current_base_vector.append("-")
                    variation = True
                    previous = '-'

            # print('abundance =', abundance/len(sequences), 'variation:', variation)
            if (abundance/len(sequences) > 0.5 and variation):
                alignment_dict[kmer] = current_base_vector

        # print("alignment dictionary python")
        # print(alignment_dict)

        vector_of_variants = []

        for i in range(0, len(paths)):
            sm = []
            for kmer in alignment_dict:
                sm.append(alignment_dict[kmer][i])
            # print(sm)
            sm.sort()
            vector_of_variants.append(sm)

        # print('to compare:', vector_of_variants)

        # C++ Code:
        # this creates three skf files which contain the k-mers + bases in binary
        for i in range(0, len(sequences)):
            # print("this is i: ", i)
            # print(sequences[i], k, paths[i], names[i])
            ska_cpp.testing_run_ska_without_reverse_complement(sequences[i], k, paths[i], names[i])

#       find skf files and read them into a vector, from binary
#       pass vector to ska_align to create alignment file from skfs
        ska_cpp.testing_run_ska_align_without_reverse_complement(paths2, names)


# test to compare each kmer that is in the alignment dictionary
        c_vector = []

        with open("/Users/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/variant_alignment.aln") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                vec = []
                for base in record.seq:
                    # print('base', base)
                    vec.append(base)
                vec.sort()
                c_vector.append(vec)

        # TODO: will not work because the ambiguity bases are not implemented in the python script
        # print('c_vector', c_vector)
        # print('vector_of_variants', vector_of_variants)
        for i in range(0, len(paths)):
            self.assertEqual(c_vector[i], vector_of_variants[i])


if __name__ == '__main__':
    unittest.main()

# TODO: test the following: reverse compliment, if Ns are skipped, ambiguous bases if they are implemented correctly
