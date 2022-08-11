#pragma once

#ifndef SKA_H
#define SKA_H

#include <vector>
#include <string>
#include "robin_hood.h"

typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;
//
//extern vec_dict kmer_dicts;

int run_ska(const std::vector<std::string>& isolate_vector, int kmerLength);

robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength);

//int run_ska_align(const std::vector<std::string>& isolate_vector, const int kmerLength);

#endif

//int run_ska() {
//    kmer_dicts(const std::vector<std::string>& isolate_vector, int kmerLength)
//}