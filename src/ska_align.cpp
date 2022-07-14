//
// Created by Johanna Helene von Wachsmann on 14/06/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include "robin_hood.h"

#include "ska.hpp"

typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;

//TODO: robin_hood::unordered_node_map better than robin_hood::unordered_map?
std::vector<robin_hood::unordered_map<std::string, char>> get_kmers(const std::vector< std::string>& isolates, const int kmer_length) {
    char split_kmer_base;
    std::string current_kmer;
    robin_hood::unordered_map<std::string, char> split_kmers;
    std::vector<robin_hood::unordered_map<std::string, char>> kmer_dicts;
    kmer_dicts.reserve(isolates.size());
    for (const auto & isolate : isolates) {
//  for (int current_isolate = 0; current_isolate < isolates.size(); current_isolate++) {
        split_kmers.clear();
        // initializing the iterators which are needed to create the kmers
//        int i1 = 0;
//      TODO: check if k-mer length is odd
        int i3 = (kmer_length - 1) / 2;
//            int i2 = i3 - 1;
        int i4 = i3 + 1;
//            int i5 = kmer_length - 1;
        int n = (kmer_length - 1) / 2;
        for (int i1 = 0; i1 <= isolate.length() - kmer_length; i1++) {
            // create kmer for dictionary
//          TODO: rolling k-mers
            current_kmer = isolate.substr(i1, n) + isolate.substr(i4, n);
            split_kmer_base = isolate[i3];
            split_kmers[current_kmer] = split_kmer_base;
            i3++; i4++; //i1++;
        }
        kmer_dicts.push_back(split_kmers);
    }

    return kmer_dicts;
}
//std::vector<char> create_bitvector_value(const std::vector<char>& bases, int dim) {
//    //TODO: turn vector of char into bitvector to be more space efficient?
////    bases = [G, C, G] and bases = [A, -, A]
//
//}

robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(std::vector<robin_hood::unordered_map<std::string, char>>& kmer_dicts) {
    robin_hood::unordered_map<std::string, std::vector<char>> all_kmers_dict;
    robin_hood::unordered_map<std::string, char> merge_dict;
    double percentage;
//    merge all keys from the single dictionaries into a vector
//    robin_hood::unordered_map<std::string, char> all_string_dict;

//    merge all the keys together into a dictionary
    for (auto & kmer_dict : kmer_dicts) {
//  for (int dictionary = 0; dictionary < kmer_dicts.size(); dictionary++) {
        merge_dict.insert(kmer_dict.begin(), kmer_dict.end());
    }
//    TODO: CHANGE vector
//  extract all the keys from dictionary to vector
    std::vector<std::string> all_keys;
    for (const auto& key: merge_dict) {
        all_keys.push_back(key.first);
    }
    double counter = 0.0;
    std::vector<char> bases;
    merge_dict.clear();
    for (auto& this_key : all_keys) {
        for (auto& kmer_dict : kmer_dicts) {
//    for (int i = 0; i < all_keys.size() ;i++) {
//        for (int j = 0; j < kmer_dicts.size(); j++) {

            if (kmer_dict.count(this_key) != 0) {
                counter++;
                bases.push_back(kmer_dict[this_key]);
            }
            else {
                bases.push_back('-');
            }
        }
        percentage = counter/kmer_dicts.size();
//      compare if current k-mer i aka all_keys[i] is in more than 90% of dictionaries of kmer_dicts
//      TODO: change percentage
        if (counter > 0 && percentage > .5) {
            std::cout << "count " << counter << " and dict size " << kmer_dicts.size() << " = " << percentage << std::endl;
            for (int b = 0; b < bases.size(); b++) {
//                create_bitvector_value(bases, kmer_dicts.size());
                all_kmers_dict[this_key] = bases;
                std::cout << "bases: " << bases[b] << " kmer: " << this_key <<std::endl;
            }
        }
        bases.clear();
        counter = 0;
    }
//printing out the contents of the dict first the key and then in the second loop the values
    for (auto const &pair: all_kmers_dict) {
        std::cout << "kmer " << pair.first << std::endl;
        for (char blub: all_kmers_dict[pair.first]){
            std::cout << "split bases " << blub << std::endl;
        }
    }
    return all_kmers_dict;
}

int run_ska(const std::vector<std::string>& isolate_vector, const int kmerLength) {

//  TODO: is this how the sequences should be inputted?
    // kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
    std::vector<robin_hood::unordered_map<std::string, char>> kmer_dicts = get_kmers(isolate_vector, kmerLength);

//  TODO: pass kmer_dicts as reference to not copy everything over all the time
//  printing vector of dictionaries
    for (int i = 0; i < kmer_dicts.size(); i++) {
        for (const auto& x: kmer_dicts[i]) {
            std::cout << "isolate dict no: " << i << ", dict size:" << kmer_dicts[i].size() << " isolate " << x.first << " "
                      << x.second << std::endl;
            // printf("isolate dict no: %d\tdict size: %lu isolate: %s %s\n", i, kmer_dicts[i].size(), x.first.c_str(), x.second);
        }
        std::cout << std::endl;
    }
    create_one_large_dictionary(kmer_dicts);
//    all_kmer_dicts = create_one_large_dictionary(kmer_dicts);
//  TODO: print kmer_dicts
//    for (auto x: kmer_dicts) {
//        std::cout << "isolate dict no: " << ", dict size:" << kmer_dicts.size() << " isolate " << x.first << " "
//                  << x.second << std::endl;
//    }
    return 0;
}

