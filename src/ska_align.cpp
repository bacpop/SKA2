//
// Created by Johanna Helene von Wachsmann on 28/07/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include "robin_hood.h"

#include "ska.hpp"
#include "ska_fasta.cpp"

//
//for (const auto & isolate : isolates)
//{
////        for (int i = 0; i < contig_count; i++){
//
//
////  for (int current_isolate = 0; current_isolate < isolates.size(); current_isolate++) {
//robin_hood::unordered_map<std::string, char> split_kmers;
//// initializing the iterators which are needed to create the kmers
////        int i1 = 0;
////      TODO: check if k-mer length is odd
//kmer_length = check_kmer_length(kmer_length);
//
//int i3 = (kmer_length - 1) / 2;
////            int i2 = i3 - 1;
//int i4 = i3 + 1;
////            int i5 = kmer_length - 1;
//int n = (kmer_length - 1) / 2;
//for (int i1 = 0; i1 <= isolate.length() - kmer_length; i1++)
//{
//// create kmer for dictionary
////          TODO: rolling k-mers
//current_kmer = isolate.substr(i1, n) + isolate.substr(i4, n);
//split_kmer_base = isolate[i3];
//split_kmers[current_kmer] = split_kmer_base;
//i3++; i4++; //i1++;
//}
////        }
//kmer_dicts.push_back(split_kmers);
//}