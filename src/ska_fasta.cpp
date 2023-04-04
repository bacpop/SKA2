//
// Created by Johanna Helene von Wachsmann on 14/06/2022.
//
#include <fstream>
#include <vector>
#include <iostream>
#include <zlib.h>
#include <chrono>
#include <sstream>
#include "robin_hood.h"
#include "ska.hpp"
#include "kseq.h"
#include <pybind11/pybind11.h>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <cereal/archives/binary.hpp>
#include "progressbar.hpp"
#include "cereal/include/cereal/archives/binary.hpp"
#include "cereal/include/cereal/types/vector.hpp"
#include "cereal/include/cereal/types/string.hpp"
#include <stdio.h>
#include <stdlib.h>
//load module openmp
KSEQ_INIT(int, read)

static const char look_up_table2[16] = {'A','C', 'G', 'T', 'M', 'R', 'W', 'S', //
                                        'Y', 'K', 'V', 'H', 'D', 'B', 'N', 'X'// }

};

//// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0b00;
static const uint64_t seedC = 0b01;
static const uint64_t seedG = 0b10;
static const uint64_t seedT = 0b11;

static const uint64_t seedM = 0b0100;
static const uint64_t seedR = 0b0101;
static const uint64_t seedW = 0b0110;
static const uint64_t seedS = 0b0111;
static const uint64_t seedY = 0b1000;
static const uint64_t seedK = 0b1001;
static const uint64_t seedV = 0b1010;
static const uint64_t seedH = 0b1011;
static const uint64_t seedD = 0b1100;
static const uint64_t seedB = 0b1101;
static const uint64_t seedN = 0b1110;
static const uint64_t seedX = 0b1111;

static const uint64_t look_up_table[256] = {
        seedX, seedT, seedX, seedG, seedA, seedA, seedX, seedC, // 0..7
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 8..15
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 16..23
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 24..31
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 32..39
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 40..47
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 48..55
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 56..63
        seedX, seedA, seedB, seedC, seedD, seedX, seedX, seedG, // 64..71
        seedH, seedX, seedX, seedK, seedX, seedM, seedN, seedX, // 72..79
        seedX, seedX, seedR, seedS, seedT, seedT, seedV, seedW, // 80..87
        seedX, seedY, seedX, seedX, seedX, seedX, seedX, seedX, // 88..95
        seedX, seedA, seedB, seedC, seedD, seedX, seedX, seedG, // 96..103
        seedH, seedX, seedX, seedK, seedX, seedM, seedN, seedX, // 104..111
        seedX, seedX, seedR, seedS, seedT, seedT, seedV, seedW, // 112..119
        seedX, seedY, seedX, seedX, seedX, seedX, seedX, seedX, // 120..127
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 128..135
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 136..143
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 144..151
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 152..159
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 160..167
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 168..175
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 176..183
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 184..191
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 192..199
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 200..207
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 208..215
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 216..223
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 224..231
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 232..239
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 240..247
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX  // 248..255
};

std::string print_kmers(int kmer_length, uint64_t kmer)
{
        std::string key = std::bitset<60>(kmer).to_string();
        std::string current_kmer("");
        for (int i = 0; i < (kmer_length-1)*2; i+=2)
        {
            std::string current_base = key.substr(i, 2);
            char base;

            if (current_base == "00")
            {
                base = 'A';
            }
            else if (current_base == "01")
            {
                base = 'C';
            }
            else if (current_base == "10")
            {
                base = 'G';
            }
            else if (current_base == "11") {
                base = 'T';
            }
//            std::cout << base;
            current_kmer.push_back(base);
        }
//        std::cout << "\n";
        return current_kmer;
}


uint64_t to_binary(std::string& current_kmer, int length) {
    // convert k-mer to bitvector
    uint64_t packed_int = 0;
    for (auto it = current_kmer.cbegin(); it != current_kmer.cend(); ++it) {
        packed_int = packed_int << 2;
        packed_int += look_up_table[*(it)];
    }
//    std::bitset<64> x(packed_int);
    return packed_int;
}

uint64_t ReverseComp64(const uint64_t mer, uint8_t kmerSize)
{
    uint64_t res = ~mer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    res = res >> (2ULL * (32 - kmerSize));
    return ((res < mer) ? res : mer);

}

uint64_t ambiguous_bases[60] = {0b0,0b0100,0b0101,0b0110,0b0100,0b0101,0b0110,0b1010,0b1011,0b1100,0b1010,0b1011,0b1100,0b1110,0b1110,
        0b0100, 0b01, 0b0111, 0b1000, 0b0100, 0b1010, 0b1011, 0b0111, 0b1000, 0b1101, 0b1010, 0b1011, 0b1110, 0b1101, 0b1110,
        0b0101, 0b0111, 0b10, 0b1001, 0b1010, 0b0101, 0b1100, 0b0111, 0b1101, 0b1001, 0b1010, 0b1110, 0b1100, 0b1101, 0b1110,
        0b0110, 0b1000, 0b1001, 0b11, 0b1011, 0b1100, 0b0110, 0b1101, 0b1000, 0b1001, 0b1110, 0b1011, 0b1100, 0b1101, 0b1110};

uint64_t no_ambiguous_bases[60] = {0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,0b00,
                                   0b01, 0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,0b01,
                                   0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10, 0b10,
                                   0b11, 0b11, 0b11, 0b11, 0b11, 0b11,0b11,0b11, 0b11,0b11,0b11,0b11,0b11,0b11, 0b11};

std::vector<int> check_for_N(std::string split, int pos)
{
    std::vector<int> results(2);
    for (int i = 0; i < split.length(); i++)
    {
        // N found
        if (split[i] == 'N')
        {
            results[0] = 1;
            // position of the N in this k-mer
            results[1] = i + 1;
            return results;
        }
    }
    return results;
}

int out_archiving(std::string path, std::vector<int64_t> test_vector)
{
    std::ofstream os(path, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(test_vector);
    return 0;
}


std::string init(int pos, int k, std::string& sequence)
{
    std::string current_kmer = sequence.substr(pos, k);
    return current_kmer;
}

int make_the_mers(int pos, std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t>& dict)
//robin_hood::unordered_map<uint64_t, uint8_t> make_the_mers(int pos, std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t>& dict)
{
    // creating k-mer mask to be used later to combine reverse compliment
    uint64_t kmer_mask = 0;
    for (int i = 0; i < k/2; i++)
    {
        kmer_mask = kmer_mask << 2;
        kmer_mask += 3;
    }

    pos = 0;
    uint64_t new_base;
    uint64_t b1;
    uint64_t b2;
    uint64_t split_base;
    uint64_t smallest_canonical;
    uint64_t mask = 3;
    // ------------------------------------
    start:
    // get first k-mer
    std::string current_kmer = init(pos, k, sequence);
    // check for Ns
    std::vector<int> check = check_for_N(current_kmer, pos);
    // an N is found
    if (check[0] == 1) {
//        std::cout << "Case 1: N found in initial k-mer." << std::endl;
        int new_pos = check[1];
        pos = pos + new_pos; //to start with new k-mer directly after the N
        if (pos < sequence.length() - k + 1) {
            goto start;
        }
    }
    // no N is found in current string k-mer
    else
    {
//        std::cout << "Case 2: no N found in initial k-mer." << std::endl;
        // convert to binary
        std::string split1 = current_kmer.substr(0, k/2);
        std::string split2 = current_kmer.substr(k/2 + 1, k/2);
        b1 = to_binary(split1, k/2);
        b2 = to_binary(split2, k/2);

        // mask the split k-mer in binary
        uint64_t binary_split_kmer = (b1 << (k/2)*2) | b2;
        // convert the split base to binary
        split_base = look_up_table[current_kmer[k/2]];

        //check for reverse compliment in binary
        smallest_canonical = ReverseComp64(binary_split_kmer, k-1);
        //
        int64_t prev_split_base = split_base;
        if (smallest_canonical < binary_split_kmer)
        {
            // if reverse compliment is smaller than current split k-mer we flip the split base
            split_base = ~split_base;
            split_base = split_base & mask;
        }
        else
        {
            // if the split k-mer is canonically smaller we continue with the split k-mer as the smallest canonical
            smallest_canonical = binary_split_kmer;
        }
        
        // check if current k-mer is already in the dictionary
        // if current split base is not in dictionary we insert it
        if (dict.find(smallest_canonical) == dict.end())
        {
            //not in dict
            uint8_t u8_split_base = split_base;
            dict[smallest_canonical] = u8_split_base;
            split_base = prev_split_base;
        }
        // if current split base is in dictionary we add the new split base based on ambiguity codes
        else
        {
            // calculate ambiguity base by calculating the offset of the vector
            int64_t u64_spilt_base = ambiguous_bases[(split_base * 15) + dict[smallest_canonical]];
//            std::cout << "ambiguity1: pos:" << pos << " " << print_kmers(k,smallest_canonical) << ": amb base: " << look_up_table2[u64_spilt_base] << " from: " << look_up_table2[split_base] << " and " << look_up_table2[dict[smallest_canonical]] << std::endl;
            uint8_t u8_split_base = u64_spilt_base;
            dict[smallest_canonical] = u8_split_base;
            split_base = prev_split_base;
        }

        // add to position to get next k-mer
        pos++;

        //getting the next k-mer from new pos
        for (int i = pos; i < sequence.length() - k + 1; i++)
        {
            new_base = look_up_table[toupper(sequence[i + k - 1])];
            if (new_base == 14)
            {
//                std::cout << "Case 3: new base is an N." << std::endl;
                //update position
                pos = i + k;
                if (pos < sequence.length() - k + 1) {
                    goto start;
                } else {
                    return 0;
                }
            }
            else
            {
//                std::cout << "Case 4: the new base is not an N." << std::endl;
                b1 = b1 << 2;
                b1 += split_base;
                b1 = b1 & kmer_mask;
                split_base = b2 >> (((k/2)*2)-2); //TODO: does this work??
                // new_base = look_up_table[sequence[k + new_base_pos]];
                b2 = b2 << 2;
                b2 += new_base;
                b2 = b2 & kmer_mask;
                uint64_t binary_split_kmer = (b1 << (k/2)*2) | b2; //bitstring
                smallest_canonical = ReverseComp64(binary_split_kmer, k-1);

                int64_t prev_split_base = split_base;
                if (smallest_canonical < binary_split_kmer)
                {
                    split_base = ~split_base;
                    split_base = split_base & mask;
                }
                else
                {
                    smallest_canonical = binary_split_kmer;
                }
                // update value if kmer already existed
                //key not present in dictionary
                if (dict.find(smallest_canonical) == dict.end())
                {
                    uint8_t u8_split_base = split_base;
                    dict[smallest_canonical] = u8_split_base;
                    split_base = prev_split_base;
                }
                else
                {
                    // calculate offset of position in vector
                    int64_t u64_spilt_base = ambiguous_bases[(split_base * 15) + dict[smallest_canonical]];
//                    std::cout << "ambiguity2: pos: " << i << ":" << print_kmers(k,smallest_canonical) << ": amb base: " << look_up_table2[u64_spilt_base] << " from: " << look_up_table2[split_base] << " and " << look_up_table2[dict[smallest_canonical]] << std::endl;
                    uint8_t u8_split_base = u64_spilt_base;
                    dict[smallest_canonical] = u8_split_base;
                    split_base = prev_split_base;
                }
            }
        }
    }
    return 0;
//    return dict;
}

void new_rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> &dict)
//robin_hood::unordered_map<uint64_t, uint8_t> new_rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> &dict)
{
    int pos = 0;
    make_the_mers(pos, sequence, k, dict);
//    return dict;
}

int check_kmer_length(int length)
{
    bool even = length % 2 == 0;
    bool max_length = length > 31;
    if (even)
    {
        std::cout << "K-mers have to be odd in length so k-mer length was changed to " << length - 1 << "!" << std::endl;
    }
    else if (max_length)
    {
        std::cout << "K-mers cannot be longer than 31!" << std::endl;
        length = 31;
    }
    return length;
}

void get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length, std::string output_directory)
//vec_dict_bits get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length, std::string output_directory)
{
    std::cout << "Ska2 fasta" << std::endl;
    // progressbar bar1(fasta_path.size());
    kmer_length = check_kmer_length(kmer_length);
    // typedef std::vector<robin_hood::unordered_map<uint64_t, uint8_t>> vec_dict_bits;
    vec_dict_bits kmer_dicts;
    kmer_dicts.resize(fasta_path.size());

    auto start = std::chrono::steady_clock::now();
    bool interrupt = false;
    std::vector<std::runtime_error> errors;

    // read in fasta files with kseq.h
//    #pragma omp parallel for // num_threads(8) // #pragma omp parallel is for paralyzing this loops
    #pragma parallel for num_threads(8)
    for (int sample_idx = 0; sample_idx < fasta_path.size(); ++sample_idx)
    {
        auto start1 = std::chrono::steady_clock::now();
//        #pragma omp critical
        // bar1.update();
        if (interrupt || PyErr_CheckSignals() != 0)
        {
            interrupt = true;
        }
        else {
            robin_hood::unordered_map<uint64_t, uint8_t> split_kmers;
            int number_of_seqs = 0;
//            #pragma omp critical
            if (!std::filesystem::exists(fasta_path[sample_idx].c_str())) {
                throw std::runtime_error("The given file does not exist!");
            }
            FILE *fp;
            kseq_t *seq;
            fp = fopen(fasta_path[sample_idx].c_str(), "r");
            seq = kseq_init(fileno(fp));
            // open current kmer file to write to string
            while (kseq_read(seq) >= 0) {
                std::string current_contig = seq->seq.s;
                ++number_of_seqs;
                new_rolling_kmer_bitvector(current_contig, kmer_length, split_kmers);
            }
            if (!interrupt) {
                // TODO: use std::move here rather than copy?
                kmer_dicts[sample_idx] = split_kmers;
                kseq_destroy(seq);
                fclose(fp);
            }
            if (split_kmers.size() == 0) {
                // #pragma omp critical
                {
                    errors.push_back(std::runtime_error("No sequence found in " + fasta_path[sample_idx]));
                    interrupt = true;
                }
            }

            auto end1 = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds = end1-start1;
            std::cout << elapsed_seconds.count() << "s\n";

            std::vector <int64_t> test_vector;

            //TODO: via index not push_back
            for (auto &it: split_kmers) {
//                std::cout << "Before stream:" << print_kmers(kmer_length, it.first) << " " << look_up_table2[it.second]
//                          << ", ";
                test_vector.push_back(it.first);
                test_vector.push_back(it.second);
            }

//            std::cout << fasta_path[sample_idx] << std::endl;
            std::string current_output_file = output_directory + names[sample_idx] + ".skf";
//            std::cout << current_output_file << std::endl;

            // binary printing as string
            // cereal printing in vector format
            out_archiving(current_output_file, test_vector);
        }
    }
//    // Check for errors
//    if (interrupt)
//    {
//        for (auto i = errors.cbegin(); i != errors.cend(); ++i)
//        {
//            std::cout << i->what() << std::endl;
//        }
//        if (errors.size())
//        {
//            throw std::runtime_error("Errors while creating k-mer dictionaries!");
//        }
//        else
//        {
//            throw pybind11::error_already_set();
//        }
//    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time total: " << elapsed_seconds.count() << "s\n";
//    std::string current_output_file = output_directory + names[0] + ".skf";
//    std::cout << "Test!" << std::endl;
//    std::cout << "Read in path: " << current_output_file << std::endl;
//    std::vector<int> dst2;
//    std::ifstream instream(current_output_file, std::ios::binary);
//    cereal::BinaryInputArchive iarchive(instream); // Create an input archive
//    iarchive(dst2);
//    std::cout << "In archive has size: " << dst2.size() << ": " << dst2[0] << std::endl;
//    for (int i = 0; i < dst2.size(); i=1+2)
//    {
//        std::cout << "dst: " << dst2[i] << ": " << print_kmers(kmer_length, dst2[i]) << std::endl;
//    }
//    return kmer_dicts;
}


std::unordered_map<uint64_t, uint8_t> change_type(robin_hood::unordered_map<uint64_t, uint8_t> rh_dict)
{
    std::unordered_map<uint64_t, uint8_t> basic_map;
    for (const auto& x: rh_dict)
    {
         basic_map[x.first] = x.second;
    }
    return basic_map;
}

robin_hood::unordered_map<uint64_t, uint8_t> change_type2(std::unordered_map<uint64_t, uint8_t> normal_dict)
{
    robin_hood::unordered_map<uint64_t, uint8_t> basic_map;
    for (const auto& x: normal_dict)
    {
        basic_map[x.first] = x.second;
    }
    return basic_map;
}

//vec_dict_bits ska_fasta(const std::vector<std::string>& isolate_vector, const std::vector<std::string>& isolate_names, int kmerLength, std::string output_directory)
int ska_fasta(const std::vector<std::string>& isolate_vector, const std::vector<std::string>& isolate_names, int kmerLength, std::string output_directory)
{
        // kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
//        vec_dict_bits kmer_dicts = get_kmers(isolate_vector, isolate_names, kmerLength, output_directory);
        get_kmers(isolate_vector, isolate_names, kmerLength, output_directory);

        // TODO: pass kmer_dicts as reference to not copy everything over all the time
//        std::cout << "n: " << kmer_dicts.size() << std::endl;
//        return kmer_dicts;
        return 0;
}


int run_ska_fasta(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength, std::string output_directory)
{
//    vec_dict_bits kmer_dicts = ska_fasta(isolate_paths, isolate_names, kmerLength, output_directory);
    ska_fasta(isolate_paths, isolate_names, kmerLength, output_directory);
    return 0;
}


std::vector<std::string> testing_run_ska(std::vector<std::string> seqs, int k, std::string kmer_or_base)
{
    // so the input are all contigs of one sequence;
    robin_hood::unordered_map<uint64_t, uint8_t> in_dict;
    robin_hood::unordered_map<uint64_t, uint8_t> out_dict;

    for (auto seq : seqs)
    {
//        std::cout << "Sequence: " << seq << std::endl;
        new_rolling_kmer_bitvector(seq, k, in_dict);
        // so the out_dict contains all k-mers of one isolate. This would be printed to the .skf file
    }
    out_dict = in_dict;
    std::vector<std::string> return_vector;
    if (kmer_or_base == "kmer")
    {
        for (auto const &pair: out_dict)
        {
            return_vector.push_back(print_kmers(k, pair.first));
        }
    }
    else if (kmer_or_base == "base")
    {
        for (auto const &pair: out_dict)
        {
            char base = look_up_table2[pair.second];
            std::string s;
            s = base;
            return_vector.push_back(s);
        }
    }
    else if (kmer_or_base == "both")
    {
        for (auto const &pair: out_dict)
        {
            char base2 = look_up_table2[pair.second];
            std::string kmer = print_kmers(k, pair.first);
            std::string s2;
            s2 = kmer + ":" + base2;
            return_vector.push_back(s2);
        }
    }
    else
    {
        std::cout << "Wrong input: valid is kmer or base!";
    }

    return return_vector;
}

const char test_ambiguity_bases(char old_base, char split_base)
{
    int index1 = look_up_table[split_base];
    int index2 = look_up_table[old_base];
    // look up in look_up_table2 to return as char and not binary int
    return look_up_table2[ambiguous_bases[(index1 * 15) + index2]];
}

std::vector<std::string> test_check_for_Ns(std::string& sequence, int k)
{
    std::vector<std::string> return_vector;
    for (int i = 0; i < sequence.length()-k+1; i++)
    {
        std::string current_kmer;
        current_kmer = sequence.substr(i, k);
        if (check_for_N(current_kmer, i)[0])  // [0] indicates if an N was found
        {
            i = check_for_N(current_kmer, i)[1] + i - 1;
        }
        else // no N was found
        {
            return_vector.push_back(current_kmer);
        }
    }
    return return_vector;
}

std::unordered_map<std::string, std::string> test_new_rolling_kmer_binary(std::string& sequence, int k)
{
//    std::cout << "sequence: " << sequence << std::endl;
    std::unordered_map<std::string, std::string> output_dict;
    robin_hood::unordered_map<uint64_t, uint8_t> integer_dict;

    new_rolling_kmer_bitvector(sequence, k, integer_dict);

    for (auto& it: integer_dict)
    {
         // TODO: only works if print bitset has k-1 size;
//         std::cout << "Key: " << print_kmers(k, it.first) << " from: " << it.first << " and value: " << look_up_table2[it.second] << " from " << it.second << std::endl;
         std::string current_kmer = print_kmers(k, it.first);
         std::string temp;
         temp = look_up_table2[it.second];
         output_dict[current_kmer] = temp;
    }
    return output_dict;
}