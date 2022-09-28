#pragma once

#ifndef SKA_H
#define SKA_H

#include <vector>
#include <string>
#include <bitset>
#include <unordered_map>
#include "robin_hood.h"

typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;
typedef std::vector<robin_hood::unordered_map<uint64_t, uint64_t>> vec_dict_bits;

int run_ska(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names,int kmerLength);

robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength);

std::vector<std::bitset<4>> create_bitvector_value(std::vector<char> bases);

uint64_t to_binary(std::string& current_kmer, int length);

robin_hood::unordered_map<uint64_t, uint64_t> rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint64_t> dict);

robin_hood::unordered_map<std::string, char> old_kmer_approach(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict);

uint64_t ReverseComp64(const uint64_t mer, uint8_t kmerSize);

std::vector<int> check_for_N(std::string split, int pos);

vec_dict_bits get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length);

//std::vector<std::unordered_map<uint64_t, uint64_t>> get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length);
#endif
