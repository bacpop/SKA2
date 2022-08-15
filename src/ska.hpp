#pragma once

#ifndef SKA_H
#define SKA_H

#include <vector>
#include <string>
#include "robin_hood.h"

typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;

int run_ska(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names,int kmerLength);

robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength);

#endif
