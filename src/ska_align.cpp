//
// Created by Johanna Helene von Wachsmann on 28/07/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include "robin_hood.h"

#include "ska.hpp"

typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;

//robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength)
robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength)
{
    robin_hood::unordered_map<std::string, std::vector<char>> all_kmers_dict;
    robin_hood::unordered_map<std::string, char> merge_dict;
    double percentage;

    // merge all the keys together into a dictionary
    for (auto & kmer_dict : kmer_dicts) {
    merge_dict.insert(kmer_dict.begin(), kmer_dict.end());
}

// the counter is used to count in how many isolates the current k-mer is present. The k-mer is discarded if a
// certain threshold isn't reached
double counter = 0.0;
std::vector<char> bases;

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
    if (counter > 0 && percentage > .5)
    {
        for (int b = 0; b < bases.size(); b++)
        {
        // create_bitvector_value(bases, kmer_dicts.size());
        all_kmers_dict[key.first] = bases;
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

return all_kmers_dict;
}

//int run_ska_align(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
////    std::cout << "ska_align.cpp" << std::endl;
////    vec_dict kmer_dicts;
//}
