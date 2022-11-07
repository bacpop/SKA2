//
// Created by Johanna Helene von Wachsmann on 28/07/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include "robin_hood.h"
#include <cereal/archives/binary.hpp>
#include "progressbar.hpp"
#include "cereal/include/cereal/archives/binary.hpp"

#include "ska.hpp"


robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> ska_align(vec_dict_bits& kmer_dicts, int kmerLength)
//robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength)
{
    robin_hood::unordered_map<std::uint64_t, std::vector<uint8_t>> alignemnt_dict;
    robin_hood::unordered_map<std::uint64_t, uint8_t> merge_dict;
    double percentage;

    // merge all the keys together into a dictionary
    for (auto & kmer_dict : kmer_dicts) {
    merge_dict.insert(kmer_dict.begin(), kmer_dict.end());
}

// the counter is used to count in how many isolates the current k-mer is present. The k-mer is discarded if a
// certain threshold isn't reached
double counter = 0.0;
std::vector<uint8_t> bases;

// iterate through the merged dictionary's keys (aka all k-mers) to check if they are also present in the
// individual dictionaries
for (auto& key: merge_dict)
{
    for (auto& kmer_dict : kmer_dicts)
    {
        if (kmer_dict.count(key.first) != 0) {
        counter++;
        bases.push_back(kmer_dict[key.first]);
        }
        else {
        // if not present in this isolate we add a placeholder
        bases.push_back('-');
        }
    }

    percentage = counter/kmer_dicts.size();
    // compare if current k-mer i aka all_keys[i] is in more than 50% of dictionaries of kmer_dicts
    // TODO: change percentage
//     bitbases;
    if (counter > 0 && percentage > .6)
    {
        for (int b = 0; b < bases.size(); b++)
        {
//            bit_bases = create_bitvector_value(bases, kmer_dicts.size());
            alignemnt_dict[key.first] = bases;
        //const int end_substring = (kmerLength-1)/2;
        //std::cout << key.first.substr(0, end_substring) << "{" << bases[b] << "}" << key.first.substr(end_substring, kmerLength) << std::endl;
        }
    }
    //    else {
    //        std::cout << "This k-mer was discarded!" << std::endl;
    //    }
    bases.clear();
    counter = 0;
}

return alignemnt_dict;
}


//int run_ska_align(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
////    std::cout << "ska_align.cpp" << std::endl;
////    vec_dict kmer_dicts;
//}
//std::vector<std::bitset<4>> create_bitvector_value(std::vector<char> bases) {
//    uint64_t packed_int = 0;
//    int base;
//    std::vector<std::bitset<4>> bitvector;
//    for (auto it = bases.cbegin(); it != bases.cend(); ++it) {
//        switch (*(it)) {
//            case 'A':
//                base = 8; // optimise out
//                break;
//            case 'C':
//                base = 4;
//                break;
//            case 'G':
//                base = 2;
//                break;
//            case 'T':
//                base = 1;
//                break;
//            case '-':
//                base = 0;
//                break;
//        }
//        std::bitset<4> x(base);
//        bitvector.push_back(x);
//    }
//
//    for (auto i: bitvector)
//        std::cout << i << ' ';
//    std::cout << "\n";
//    return bitvector;
//}
vec_dict_bits ska_merge() {
//TODO!
}

int run_ska_align(const std::vector<std::string>& skf_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
    std::string line;
    int blub = 0;
//    TODO: parallel
    for (std::string file : skf_paths){
        std::ifstream instream(file, std::ios::binary);
        cereal::BinaryInputArchive iarchive(instream); // Create an input archive
        MyBase cereal_split_kmers;
//        cereal_split_kmers.cereal_dict;
        iarchive(cereal_split_kmers.cereal_dict);
//        if (instream.is_open())
//        {
//            while (getline (instream,line)) {
//                blub++;
//                if (blub % 1000 == 0){
//                    std::cout << blub << std::endl;
//                }
//            }
//        }
//        robin_hood::unordered_map<uint64_t, uint8_t> test_dict = cereal_split_kmers;
//        std::cout << test_dict.size() << std::endl;
    }


//    std::ofstream ins("/Users/wachsmannj/Documents/test_SKA2/integer_approach/testing/" +names[sample_idx] + ".skf", std::ios::binary);
//    // creating archive
//    cereal::BinaryOutputArchive oarchive(os);
//    MyBase cereal_split_kmers;
//    cereal_split_kmers.cereal_dict = split_kmers;
//    oarchive(cereal_split_kmers);

//    cereal::BinaryInputArchive iarchive(ins); // Create an input archive
//    MyBase cereal_split_kmers;
//    cereal_split_kmers.cereal_dict = split_kmers;
//    iarchive(cereal_split_kmers);

    std::cout << "blub" << std::endl;
    //    do something to create k-mer dict, this is just reading in all the skf and putting them into a vector
//    vec_dict_bits kmer_dicts = ska_merge();
//    then the vector can be put into the actual alignment creating function
//    robin_hood::unordered_map<std::uint64_t, std::vector<uint8_t>> alignment = ska_align(kmer_dicts, 31);

}