import unittest
import ska_cpp
import math
from Bio import SeqIO
from ska.__main__ import read_in_files


def reverse_compliment(seq):
    comp = []
    for base in seq:
        if base == "A":
            comp.append("T")
        elif base == "G":
            comp.append("C")
        elif base == "T":
            comp.append("A")
        elif base == "C":
            comp.append("G")

    # reverse the sequence
    comp_rev = comp[::-1]

    # convert list to string
    comp_rev = "".join(comp_rev)
    return comp_rev


class MyTestCase(unittest.TestCase):

    def test_ska_align(self):
        k = 7

        rfile = "/home/johanna/test_SKA2/integer_approach/MA_results/rfile.txt"
        output_dir = "/home/johanna/test_SKA2/integer_approach/MA_results/"
        # rfile = "/home/johanna/test_SKA2/test_run_ska/rfile_isolates.txt"
        # output_dir = "/home/johanna/test_SKA2/test_run_ska/"
        cluster_name = "cluster28"
        paths, names = [], []
        with open(rfile, newline='\n') as file:
            lines = file.readlines()
            for line in lines:
                paths.append(line.split("\t")[1].rstrip())
                names.append(line.split("\t")[0].rstrip())

        # run C++ code
        ska_cpp.run_ska_fasta(paths, names, k, output_dir)

        skf_paths = []
        for sample in names:
            skf_files = output_dir + sample + ".skf"
            skf_paths.append(skf_files)
        # ska_cpp.testing_run_ska_align_without_reverse_complement(skf_paths, names, k)
        ska_cpp.run_ska_align(skf_paths, names, k, output_dir, cluster_name)
        ska_seqs = []
        aln_path = output_dir + cluster_name + ".aln"

        # python code to get alignment with reverse compliment
        # read in fasta file sequences and get them in the following format:
        # sequences [['AACAA', 'CCCCCCC', 'GGGGATATATATAG'], ['AACAA', 'CGCCC', 'GGGAG'], ['AACAA', 'CCTCC', 'GGGATATATAGGGGG']]

        with open(aln_path, 'r') as fasta_sequences:
            for fasta in SeqIO.parse(fasta_sequences, "fasta"):
                sequence = str(fasta.seq)
                ska_seqs.append(sequence)

        # pre-processing
        sequences = []
        for path in paths:
            current_sequence = []
            with open(path, 'r') as fasta_sequences:
                for fasta in SeqIO.parse(fasta_sequences, "fasta"):
                    current_contig = str(fasta.seq)
                    current_sequence.append(current_contig)
                sequences.append(current_sequence)
        # print("sequence", sequences)

        # actual code starting here
        position = 0
        vector_of_dicts = []
        c_dict = []
        for i in range(0, len(sequences)):
            current_dict = {}
            # print("sequences",  sequences[i])
            print("path:", paths[i])
            current_dict = {}
            c_mers = ska_cpp.testing_run_ska(sequences[i], k, "both")
            for j in range(0, len(sequences[i])):
                current_sequence = sequences[i][j]
                current_dict = dict_rolling_kmers(current_sequence, k, current_dict)

            vector_of_dicts.append(current_dict)

        # print('vector of dicts', vector_of_dicts)
        # create unique_merge_vector to have all unique kmers to iterate through
        merge_vector = []
        for dic in vector_of_dicts:
            merge_vector.append(list(dic.keys()))
        flat_merge_vector = [item for blub in merge_vector for item in blub]
        myset = set(flat_merge_vector)
        unique_merge_vector = list(myset)
        unique_merge_vector.sort()
        # print(unique_merge_vector)

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
                    if previous != 'J' and previous != isolate_dict[kmer]:
                        variation = True
                    previous = isolate_dict[kmer]
                else:
                    current_base_vector.append("-")
                    variation = True
                    previous = '-'
            # print(abundance, len(sequences))
            if abundance / len(sequences) > 0.9 and variation:
                alignment_dict[kmer] = current_base_vector
        vector_of_variants = []

        for i in range(0, len(paths)):
            sm = []
            for kmer in alignment_dict:
                sm.append(alignment_dict[kmer][i])
            sm.sort()
            vector_of_variants.append(sm)

        c_vector = []

        for i in range(0, len(paths)):
            sm = []
            for j in range(0, len(ska_seqs[0])):
                # print(ska_seqs[i][j])
                sm.append(ska_seqs[i][j])
            sm.sort()
            c_vector.append(sm)

        for i in range(0, len(ska_seqs)):
            print(sorted(c_vector[i]))
            print(sorted(vector_of_variants[i]))
            self.assertEqual(''.join(sorted(c_vector[i])), ''.join(sorted(vector_of_variants[i])))  # add assertion here


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


def rolling_kmers(sequence, k):
    pos = 0
    current_kmer_dict = {}
    current_kmer_dict = make_the_mers_ambiguity(pos, sequence, k, current_kmer_dict)
    return current_kmer_dict


def dict_rolling_kmers(sequence, k, current_kmer_dict):
    pos = 0
    current_kmer_dict = make_the_mers_ambiguity(pos, sequence, k, current_kmer_dict)
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

# TODO: test the following: reverse compliment, if Ns are skipped, ambiguous bases if they are implemented correctly

