from Bio.Seq import Seq
import unittest
import ska_cpp


class MyTestCase(unittest.TestCase):
    def test_ska_fasta(self):

        # Test 1: test binary conversion
        k = 7
        s = "ACTGCTA"
        bitstring = ska_cpp.to_binary(s, k)
        print("bitstring", bitstring)
        self.assertEqual(1948, bitstring)

        #         test reverse complement
        reverse = ska_cpp.ReverseComp64(bitstring, k)
        self.assertEqual(12875, reverse)
        # 00110010
        print("Reverse", reverse)

        #         test check_for_N

        check = ska_cpp.check_for_N(s, 0)
        self.assertEqual([0, 0], check)
        print("Check", check)

        # Test 2:
        k = 7
        # s = "ACTGATCCNCTT"
        s = "ACTG"
        #         test check_for_N
        check = ska_cpp.check_for_N(s, 0)
        # self.assertEqual([1, 9], check)
        print("Check", check)

        # testing ska fasta

        s = "CTAACCACCGTGTATTCGTTATGGCACCAGGGAGTTTAAGCCGAGTCAATGGAGCTCGCAATACAGAGTTTACCGCATCTTGCCCTAACTGACAAACTGT"
        k = 3
        ska_fasta = ska_cpp.testing_run_ska(s, k)
        print(ska_fasta)


        python_ska_fasta = make_kmers(s,k)

        python_ska_fasta = set(python_ska_fasta)
        s_sorted = dict.fromkeys(sorted(python_ska_fasta))
        python_ska_fasta = (list(s_sorted))
        print(python_ska_fasta)
        # python_ska_fasta = python_ska_fasta.sort()
        # print(python_ska_fasta)


        self.assertEqual(python_ska_fasta, ska_fasta)

#

def make_kmers(s,k):
    sub_len = int(((k-1)/2))
    python_ska_fasta = []
    for i in range(0, len(s)-k+1):
        if s[i] == 'N' or s[i] == 'n':
            i += k
        current_kmer = s[i: i+sub_len] + s[i+sub_len+1: i+sub_len+1+sub_len]

        reverse_kmer = Seq(current_kmer).reverse_complement()
        # accounting for reverse compliment
        if (reverse_kmer < current_kmer):
            python_ska_fasta.append(reverse_kmer)
            print(i, current_kmer)
        else:
            python_ska_fasta.append(current_kmer)
            print(i, current_kmer)
    return python_ska_fasta


if __name__ == '__main__':
    unittest.main()
