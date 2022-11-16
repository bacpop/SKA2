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
#include "cereal/include/cereal/types/vector.hpp"
#include "cereal/include/cereal/types/string.hpp"
#include "eigen/Eigen/Eigen"
#include "ska.hpp"


robin_hood::unordered_map<std::uint64_t , std::vector<uint64_t>> ska_align(std::vector<robin_hood::unordered_map<uint64_t, uint64_t>>& kmer_dicts, int kmerLength)
//robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength)
{
    int c = 0;
    robin_hood::unordered_map<int, int> base_count;
    robin_hood::unordered_map<std::uint64_t, std::vector<uint64_t>> alignemnt_dict;
    robin_hood::unordered_map<std::uint64_t, uint64_t> merge_dict;
    double percentage;

    // merge all the keys together into a dictionary
    for (auto & kmer_dict : kmer_dicts) {
    merge_dict.insert(kmer_dict.begin(), kmer_dict.end());
}

    // the counter is used to count in how many isolates the current k-mer is present. The k-mer is discarded if a
    // certain threshold isn't reached
    double counter = 0.0;
    std::vector<uint64_t> bases;

    // iterate through the merged dictionary's keys (aka all k-mers) to check if they are also present in the
    // individual dictionaries
    std::cout << "merge dict size: " << merge_dict.size() << std::endl;

    for (auto& key: merge_dict)
    {
        for (auto& kmer_dict : kmer_dicts)
        {
            if (kmer_dict.count(key.first) != 0)
            {
                counter++;
                bases.push_back(kmer_dict[key.first]);
            }
            else
            {
                // if not present in this isolate we add a placeholder
                bases.push_back('-');
            }
        }

    percentage = counter/kmer_dicts.size();
    // compare if current k-mer i aka all_keys[i] is in more than 50% of dictionaries of kmer_dicts
    // TODO: change percentage
//     bitbases;


    if (counter > 0 && percentage > .9)
    {
        for (int b = 0; b < bases.size(); b++)
        {
            alignemnt_dict[key.first] = bases;

        }
    }

//    std::cout << "Bases count: " << base_count.size() << std::endl;
    //    else {
    //        std::cout << "This k-mer was discarded!" << std::endl;
    //    }
    bases.clear();
    counter = 0;
}



return alignemnt_dict;
}


vec_dict_bits ska_merge() {
//TODO???
}

int run_ska_align(const std::vector<std::string>& skf_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
//    TODO: parallel
    std::vector<robin_hood::unordered_map<uint64_t, uint64_t>> all_kmers;

//    robin_hood::unordered_map<uint64_t, uint8_t> current_isolate_dict;
    for (std::string file : skf_paths) {
//        std::cout << "blub: " << file << std::endl;
        robin_hood::unordered_map<uint64_t, uint64_t> current_isolate_dict;
        std::vector<int> dst;
        std::ifstream instream(file, std::ios::binary);
        cereal::BinaryInputArchive iarchive(instream); // Create an input archive

        iarchive(dst);
//        std::cout << "vecotr size: " << dst.size()/2 << std::endl;

        for (int i = 0; i < dst.size(); i += 2) {
            current_isolate_dict[i] = i+1;
        }
        all_kmers.push_back(current_isolate_dict);
    }

    std::cout << all_kmers[0].size() << std::endl;

    robin_hood::unordered_map<std::uint64_t , std::vector<uint64_t>> variant_alignemnt;
    variant_alignemnt = ska_align(all_kmers, 31);

    //  TODO: make vectors for each file
    //  TODO: check if isolate names and vactor order are correct
    //  TODO: cerealize the output to file as chars

    std::cout << "printing alignments to file!" << std::endl;
    std::ofstream myfile;
    myfile.open ("/Users/wachsmannj/Documents/test_SKA2/integer_approach/small_test/test_cluster.aln");
    int i = 0;
    int c = 0;

    Eigen::MatrixXf mymatrix(all_kmers[0].size(), skf_paths.size());
    Eigen::Vector
    for (int i = 0; i < isolate_names.size(); i ++)
    {

        for (auto& it: variant_alignemnt)
        {
//            mymatrix << it.second;
            myfile << it.second[i] << "\n";
            c++;
        }

    }
    std::cout << "c: " << c<< std::endl;
    myfile.close();





}