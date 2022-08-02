//
// Created by Johanna Helene von Wachsmann on 14/06/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include "robin_hood.h"

#include "ska.hpp"
#include "kseq.h"
//KSEQ_INIT(gzFile, gzread)

typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;

bool check_kmer_length(int length) {
    bool odd = length % 2 == 1;
    if (odd) {
        std::cout << "K-mers have to be odd in length so k-mer length was changed to " << length - 1 << "!"
                  << std::endl;
    }
    return odd;
}

//TODO: robin_hood::unordered_node_map better than robin_hood::unordered_map?
vec_dict get_kmers(const std::vector< std::string>& isolates, int kmer_length, const std::vector< int>& contig_count)
{
    char split_kmer_base;
    std::string current_kmer;
    vec_dict kmer_dicts;
    kmer_dicts.reserve(isolates.size());


    int current_pos = 0;
    int isolate_length_pos = 0;

    while (current_pos < isolates.size()){
        for (int i = 0; i < contig_count.size(); i++){
            std::string current_isolate = isolates[current_pos + 1];
            std::cout << current_isolate << std::endl;




        }
        current_pos = current_pos + current_pos + 1;
        isolate_length_pos ++;
    }
//
//    for (const auto & isolate : isolates)
//    {
////        for (int i = 0; i < contig_count; i++){
//
//
////  for (int current_isolate = 0; current_isolate < isolates.size(); current_isolate++) {
//        robin_hood::unordered_map<std::string, char> split_kmers;
//        // initializing the iterators which are needed to create the kmers
////        int i1 = 0;
////      TODO: check if k-mer length is odd
//        kmer_length = check_kmer_length(kmer_length);
//
//        int i3 = (kmer_length - 1) / 2;
////            int i2 = i3 - 1;
//        int i4 = i3 + 1;
////            int i5 = kmer_length - 1;
//        int n = (kmer_length - 1) / 2;
//        for (int i1 = 0; i1 <= isolate.length() - kmer_length; i1++)
//        {
//            // create kmer for dictionary
////          TODO: rolling k-mers
//            current_kmer = isolate.substr(i1, n) + isolate.substr(i4, n);
//            split_kmer_base = isolate[i3];
//            split_kmers[current_kmer] = split_kmer_base;
//            i3++; i4++; //i1++;
//        }
////        }
//        kmer_dicts.push_back(split_kmers);
//    }

    return kmer_dicts;
}


//std::vector<char> create_bitvector_value(const std::vector<char>& bases, int dim) {
//TODO: turn vector of char into bitvector to be more space efficient?
////    bases = [G, C, G] and bases = [A, -, A]
//}



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

//    std::cout << "---- The current k-mer is: " << key.first << std::endl;
    percentage = counter/kmer_dicts.size();
    // compare if current k-mer i aka all_keys[i] is in more than 50% of dictionaries of kmer_dicts
    // TODO: change percentage
    if (counter > 0 && percentage > .5)
    {
        for (int b = 0; b < bases.size(); b++)
        {
            // create_bitvector_value(bases, kmer_dicts.size());
            all_kmers_dict[key.first] = bases;
            const int end_substring = (kmerLength-1)/2;
            std::cout << key.first.substr(0, end_substring) << "{" << bases[b] << "}" << key.first.substr(end_substring, kmerLength) << std::endl;
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



int run_ska(const std::vector<std::string>& isolate_vector, const std::vector<int>& contig_count, int kmerLength) {

//    const std::vector< int> contig_count;
    // kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
    vec_dict kmer_dicts = get_kmers(isolate_vector, kmerLength, contig_count);
//    vec_dict get_kmers(const std::vector< std::string>& isolates, int kmer_length, const std::vector< int>& contigCount)
// TODO: pass kmer_dicts as reference to not copy everything over all the time
//// printing vector of dictionaries
//    for (int i = 0; i < kmer_dicts.size(); i++)
//    {
//        for (const auto& x: kmer_dicts[i])
//        {
//            std::cout << "isolate dict no: " << i << ", dict size:" << kmer_dicts[i].size() << " isolate " << x.first << " "
//                      << x.second << std::endl;
//            // printf("isolate dict no: %d\tdict size: %lu isolate: %s %s\n", i, kmer_dicts[i].size(), x.first.c_str(), x.second);
//        }
//        std::cout << std::endl;
//    }

    create_one_large_dictionary(kmer_dicts, kmerLength);
// all_kmer_dicts = create_one_large_dictionary(kmer_dicts);
//  TODO: print kmer_dicts
//    for (auto x: kmer_dicts) {
//        std::cout << "isolate dict no: " << ", dict size:" << kmer_dicts.size() << " isolate " << x.first << " "
//                  << x.second << std::endl;
//    }
    return 0;
}
