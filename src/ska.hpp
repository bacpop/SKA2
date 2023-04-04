#pragma once

//#ifndef SKA_H
//#define SKA_H

#include <vector>
#include <string>
#include <bitset>
#include <unordered_map>
#include "robin_hood.h"
#include "robin_hood_cereal.h"


typedef std::vector<robin_hood::unordered_map<std::string, char>> vec_dict;
typedef std::vector<robin_hood::unordered_map<uint64_t, uint8_t>> vec_dict_bits;

int run_ska_fasta(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names,int kmerLength, std::string output_directory);
//int run_ska_fasta(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names,int kmerLength, std::string output_directory);

void run_ska_align(const std::vector<std::string>& skf_paths, const std::vector<std::string>& isolate_names, int kmerLength, std::string output_directory, std::string cluster_name);


//robin_hood::unordered_map<std::string, std::vector<char>> create_one_large_dictionary(vec_dict& kmer_dicts, int kmerLength);

//std::vector<std::bitset<4>> create_bitvector_value(std::vector<char> bases);

uint64_t to_binary(std::string& current_kmer, int length);

//robin_hood::unordered_map<uint64_t, uint8_t> rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> dict);

//robin_hood::unordered_map<std::string, char> old_kmer_approach(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict);

uint64_t ReverseComp64(const uint64_t mer, uint8_t kmerSize);

std::vector<int> check_for_N(std::string split, int pos);

//vec_dict_bits get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length, std::string output_directory);
void get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length, std::string output_directory);

std::unordered_map<uint64_t, uint8_t> change_type(robin_hood::unordered_map<uint64_t, uint8_t> rh_dict);

robin_hood::unordered_map<uint64_t, uint8_t> change_type2(std::unordered_map<uint64_t, uint8_t> normal_dict);
// cereal needs to know which data members to serialize in your classes. Let it know by implementing a serialize method in your class

std::vector<std::string> testing_run_ska(std::vector<std::string>  seqs, int k, std::string kmer_or_base);

//void testing_run_ska_without_reverse_complement(std::vector<std::string> seqs, int k, std::string path, std::string file_name);
//void testing_run_ska_align_without_reverse_complement(std::vector<std::string> skf_paths, std::vector<std::string> isolate_names, int k);
//std::unordered_map<uint64_t, uint8_t> testing_rolling_kmer_bitvector(std::string& sequence, int k, std::unordered_map<uint64_t, uint8_t> dict);


//TESTING FUNCTIONS:
const char test_ambiguity_bases(char old_base, char split_base);
std::vector<std::string> test_check_for_Ns(std::string& sequence, int k);
//robin_hood::unordered_map<uint64_t, uint8_t> new_rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> &dict);
void new_rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> &dict);
std::unordered_map<std::string, std::string> test_new_rolling_kmer_binary(std::string& sequence, int k);
std::string print_kmers2(int kmer_length, uint64_t kmer);

std::vector<std::string>compare_kmers(std::string file, int k);
//struct MyBase
//{
//    robin_hood::unordered_map<uint64_t, uint8_t> cereal_dict;
//
//    template <class Archive>
//    void serialize( Archive & ar )
//    {
//        ar(cereal_dict);
//    }
//};
//
//struct MyDerived : MyBase
//{
//    template <class Archive>
//    void load( Archive & )
//    { }
//
//    template <class Archive>
//    void save( Archive & ) const
//    { }
//};
//

//#endif
