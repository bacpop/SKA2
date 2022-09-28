import unittest
import ska_cpp


class MyTestCase(unittest.TestCase):
    def test_ska_fasta(self):
        # Test 1:
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

        unordered_map = ska_cpp.get_kmers(["tests/test_string.fasta"], ["test1"], 3)
        print(unordered_map)
#          whole algorithm
#         ska_cpp.rolling_kmer_extract(s, k, )





if __name__ == '__main__':
    unittest.main()
