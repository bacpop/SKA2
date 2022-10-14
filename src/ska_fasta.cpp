//
// Created by Johanna Helene von Wachsmann on 14/06/2022.
//
//#include "boost/filesystem.hpp"
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

//using namespace boost::filesystem;
KSEQ_INIT(gzFile, gzread)

//// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0b00;
static const uint64_t seedC = 0b01;
static const uint64_t seedG = 0b10;
static const uint64_t seedT = 0b11;
//static const uint64_t seedN = 0b100; // this won't be used

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

uint64_t to_binary(std::string& current_kmer, int length) {
    // convert k-mer to bitvector
    uint64_t packed_int = 0;
    for (auto it = current_kmer.cbegin(); it != current_kmer.cend(); ++it) {
        packed_int = packed_int << 2;
        packed_int += look_up_table[*(it)];
    }
    std::bitset<64> x(packed_int);
    return packed_int;
}

//// from Simon Hariss's original code
void ascii_bitstring(std::string & mybits){
    int myremainder=::fmod(int(mybits.length()),6);
    if (myremainder>0){
        for (int i = 0; i<(6-myremainder); ++i){
            mybits.push_back('0');
        }
    }
    for (std::string::size_type i = 0; i<mybits.length(); i+=6){
        assert((i+5)<mybits.length());
        mybits[i/6] = ((int(mybits[i])-'0')+((int(mybits[i+1])-'0')*2)+((int(mybits[i+2])-'0')*4)+((int(mybits[i+3])-'0')*8)+((int(mybits[i+4])-'0')*16)+((int(mybits[i+5])-'0')*32))+33;
        }
    mybits.erase(mybits.length()/6);
}

uint64_t ReverseComp64(const uint64_t mer, uint8_t kmerSize)
{
    uint64_t res = ~mer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    return(res >> (2ULL * (32 - kmerSize)));
//    res = res >> (2ULL * (32 - kmerSize));
//
//    return ((res < mer) ? res : mer);
// TODO: fix and return smallest one; and fix code further down
}

uint64_t ambiguous_bases_old[60] = {0b0, 0b0100, 0b0101, 0b0110,     //{{'A', 'M', 'R', 'W',},
                                   0b0100, 0b01, 0b0111, 0b1000,    //{'M', 'C', 'S', 'Y',},
                                   0b0101, 0b0111, 0b10, 0b1001,    //{'R', 'S', 'G', 'K',},
                                   0b0110, 0b1000, 0b1001, 0b11,    //{'W', 'Y', 'K', 'T',},
                                   0b0100, 0b0100, 0b1010, 0b1011,  //{'M', 'M', 'V', 'H',},
                                   0b0101, 0b1010, 0b0101, 0b1100,  //{'R', 'V', 'R', 'D',},
                                   0b0110, 0b1011, 0b1100, 0b0110,  //{'W', 'H', 'D', 'W',},
                                   0b1010, 0b0111, 0b0111, 0b1101,  //{'V', 'S', 'S', 'B',},
                                   0b1011, 0b1000, 0b1101, 0b1000,  //{'H', 'Y', 'B', 'Y',},
                                   0b1100, 0b1101, 0b1001, 0b1001,  //{'D', 'B', 'K', 'K',},
                                   0b1010, 0b1010, 0b1010, 0b1110, //{'V', 'V', 'V', 'N',},
                                   0b1011, 0b1011, 0b1110, 0b1011, //{'H', 'H', 'N', 'H',},
                                   0b1100, 0b1110, 0b1100, 0b1100, //{'D', 'N', 'D', 'D',},
                                   0b1110, 0b1101, 0b1101, 0b1101, //{'N', 'B', 'B', 'B',},
                                   0b1110, 0b1110, 0b1110, 0b1110};//{'N', 'N', 'N', 'N'}};

uint64_t ambiguous_bases[60] = {0b0,0b0100,0b0101,0b0110,0b0100,0b0101,0b0110,0b1010,0b1011,0b1100,0b1010,0b1011,0b1100,0b1110,0b1110,
        0b0100, 0b01, 0b0111, 0b1000, 0b0100, 0b1010, 0b1011, 0b0111, 0b1000, 0b1101, 0b1010, 0b1011, 0b1110, 0b1101, 0b1110,
        0b0101, 0b0111, 0b10, 0b1001, 0b1010, 0b0101, 0b1100, 0b0111, 0b1101, 0b1001, 0b1010, 0b1110, 0b1100, 0b1101, 0b1110,
        0b0110, 0b1000, 0b1001, 0b11, 0b1011, 0b1100, 0b0110, 0b1101, 0b1000, 0b1001, 0b1110, 0b1011, 0b1100, 0b1101, 0b1110};


robin_hood::unordered_map<std::string, char> old_kmer_approach(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict) {
    // initializing the iterators which are needed to create the kmers
//              int i1 = 0;
    int i3 = (k - 1) / 2;
//              int i2 = i3 - 1;
    int i4 = i3 + 1;
//              int i5 = k - 1;
    int n = (k - 1) / 2;
    char split_kmer_base;
    std::string current_kmer;
    for (int i1 = 0; i1 < sequence.length() - k; i1++) {
//          TODO: rolling k-mers
        current_kmer = sequence.substr(i1, n) + sequence.substr(i4, n);
        split_kmer_base = sequence[i3];
        dict[current_kmer] = split_kmer_base;
//        std::cout << current_kmer << ": " << split_kmer_base << std::endl;
        i3++;
        i4++;
    }
    return dict;
}

//robin_hood::unordered_map<std::string, char> rolling_kmer_extract(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict) {
//    std::string temp;
//    std::string current_kmer;
//    std::string split_kmer;
//    char base;
//    for (int i = k-1; i < sequence.length();  i++) {
//        if (temp != "") {
//            current_kmer = temp.substr(1, k-1) + sequence[i];
//        }
//        else {
//            current_kmer = sequence.substr(0, k);
//        }
////        base = extract_middlebase(current_kmer);
//        base = current_kmer[(k/2)];
//        temp = current_kmer;
//        split_kmer = current_kmer.erase((k/2), 1);
//        dict[split_kmer] = base;
////        std::cout << current_kmer << ": " << base << std::endl;
//    }
//    return dict;
//}

robin_hood::unordered_map<std::string, char> rolling_kmer_iterators(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict) {
    std::string temp;
    std::string current_kmer;
    int i3 = (k - 1) / 2;
    int i2 = i3 - 1;
    int i4 = i3 + 1;
    int i5 = k - 1;
    int n = (k - 1) / 2;
    for (int i1 = 0; i1 <= sequence.length() - k; i1++) {
        char split_kmer_base;
        if (temp != "") {
            current_kmer = temp.substr(i1, k-1) + sequence[i2] + temp.substr(i4, k-1) + sequence[i5];
        }
        else {
            current_kmer = sequence.substr(i1, n) + sequence.substr(i4, n);
        }
        split_kmer_base = sequence[i3];
        dict[current_kmer] = split_kmer_base;
        i3++;
        i4++;
    }
    return dict;
}

std::vector<int> check_for_N(std::string split, int pos){
    std::vector<int> results(2);
    for (int i = 0; i < split.length(); i++) {
        // N found
        if (split[i] == 'N'){
            results[0] = 1;
            // position of the N
            results[1] = i + 1 + pos;
            return results;
        }
        // N not found
        else {
            results[0] = 0;
            results[1] = pos;
        }
    }
    return results;
}


robin_hood::unordered_map<uint64_t, uint64_t> rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint64_t> dict) {
//std::unordered_map<uint64_t, uint64_t> rolling_kmer_bitvector(std::string& sequence, int k, std::unordered_map<uint64_t, uint64_t> dict) {
    uint64_t kmer_mask = 0;
    for (int i = 0; i < k/2; i++) {
        kmer_mask = kmer_mask << 2;
        kmer_mask += 3;

    }
    int pos = 0;
    bool new_base_N = true;
    uint64_t new_base;
    uint8_t sub_kmer_length = k/2;
    std::string current_kmer = sequence.substr(pos, k);
    std::vector<int> check;
    uint64_t b1;
    uint64_t b2;
    uint64_t m;
    uint64_t smallest_canonical;
    uint64_t mask = 3;

    while (new_base_N && pos < sequence.length()) {
        check = check_for_N(current_kmer, pos);
        // no mpre Ns
        if (check[0] == 0) {
            pos = check[1]; //9
            new_base = look_up_table[sequence[k + pos]] ;
            new_base_N = false;
            // convert to binary
            std::string split1 = current_kmer.substr(0, k/2);
            std::string split2 = current_kmer.substr(k/2 + 1, k/2);

            b1 = to_binary(split1, k/2);
            b2 = to_binary(split2, k/2);

            uint64_t bint = (b1 << ((k/2)+1)) | b2;
            m = look_up_table[sequence[sub_kmer_length]];
            smallest_canonical = ReverseComp64(bint, k-1);
            int64_t value = m;

            if (smallest_canonical < bint) {
                value = ~m;
                value = value & mask;
            }
            else {
                smallest_canonical = bint;
            }
            // update value if kmer already existed
            //key not present in dictionary
            uint8_t mid_val = static_cast<uint8_t>(value);
            if (dict.find(smallest_canonical) == dict.end()) {
                dict[smallest_canonical] = value;
            }
            else {
                // calculate offset of position in vector
                int64_t value = ambiguous_bases[(m * 15) + dict[smallest_canonical]];
                dict[smallest_canonical] = value;
            }
        }
            // more Ns
        else {
            pos = check[1];
            current_kmer = sequence.substr(pos, k);
            new_base_N = true;
        }
    }
//TODO: does N still = 4 here??
    new_base = look_up_table[sequence[k + pos]]; //should be either 00,01,10,11
    for (int new_base_pos = pos; new_base_pos < sequence.size() - k + 1; ++new_base_pos){
        new_base = look_up_table[sequence[k + new_base_pos]];

        if(new_base == 14) {
            new_base_N = true;
            new_base_pos = new_base_pos + k + 1;
            current_kmer = sequence.substr(new_base_pos, k);
            while (new_base_N && new_base_pos < sequence.length()) {
                check = check_for_N(current_kmer, new_base_pos);
                // no mpre Ns
                if (check[0] == 0) {
                    new_base_pos = check[1]-1; //because of loop?
                    new_base_N = false;
                }
                    // more Ns
                else {
                    new_base_pos = check[1];
                    current_kmer = sequence.substr(new_base_pos, k);
                    new_base_N = true;
                }
            }
        }
        else {
            b1 = b1 << 2;
            b1 += m;
            b1 = b1 & kmer_mask;
            m = b2 >> (((k/2)*2)-2); //does this work??
            new_base = look_up_table[sequence[k + new_base_pos]];
            b2 = b2 << 2;
            b2 += new_base;
            b2 = b2 & kmer_mask;
            uint64_t bitstring = (b1 << (k/2+1)) | b2;
            smallest_canonical = ReverseComp64(bitstring, k/2+1);
            int64_t value = m;
            if (smallest_canonical < bitstring) {
                value = ~m;
                value = value & mask;
            }
            else {
                smallest_canonical = bitstring;
            }
            // update value if kmer already existed
            //key not present in dictionary
            if (dict.find(smallest_canonical) == dict.end()) {
                dict[smallest_canonical] = value;
            }
            else {
                // calculate offset of position in vector
                int64_t value1 = ambiguous_bases[(value * 15) + dict[smallest_canonical]];
                dict[smallest_canonical] = value1;
            }
        }
    }
    return dict;
}


int check_kmer_length(int length) {
    bool even = length % 2 == 0;
    bool max_length = length > 31;
    if (even) {
        std::cout << "K-mers have to be odd in length so k-mer length was changed to " << length - 1 << "!" << std::endl;
    }
    else if (max_length) {
        std::cout << "K-mers cannot be longer than 31!" << std::endl;
        length = 31;
    }
    return length;
}

//TODO: robin_hood::unordered_node_map better than robin_hood::unordered_map?
vec_dict_bits get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length)
//std::vector<std::unordered_map<uint64_t, uint64_t>> get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length)
{
//    progressbar bar1(fasta_path.size());
//    std::cout << "Creating k-mer dictionary: " << std::endl;
    kmer_length = check_kmer_length(kmer_length);
    vec_dict_bits kmer_dicts;
//    std::vector<std::unordered_map<uint64_t, uint64_t>> kmer_dicts;
    kmer_dicts.resize(fasta_path.size());
// read in fasta files with kseq.h
    auto start = std::chrono::steady_clock::now();
    bool interrupt = false;
    std::vector<std::runtime_error> errors;


#pragma omp parallel for // num_threads(2) // #pragma omp parallel is for paralyzing this loops
    for (int sample_idx = 0; sample_idx < fasta_path.size(); ++sample_idx) {
        #pragma omp critical
//        bar1.update();
        if (interrupt || PyErr_CheckSignals() != 0) {
            interrupt = true;
        } else {
            //#pragma omp critical
            robin_hood::unordered_map<uint64_t, uint64_t> split_kmers;
//            std::unordered_map<uint64_t, uint64_t> split_kmers;
            gzFile fp;
            kseq_t *seq;
            int number_of_seqs = 0;
//          TODO: change filesystem to boost
//          boost::filesystem::exists(fasta_path[sample_idx]);
            std::filesystem::path p(fasta_path[sample_idx]);
            if (!std::filesystem::exists(p)) {
                throw std::runtime_error("The given file does not exist!");
            }
            fp = gzopen(fasta_path[sample_idx].c_str(), "r");
            seq = kseq_init(fp);
            //open current kmer file to write to string
            while ((kseq_read(seq)) >= 0) {
                std::string current_contig = seq->seq.s;
                ++number_of_seqs;
//              TODO: is contig dict updated or is dict only for contigs???

//                split_kmers = rolling_kmer_iterators(current_contig, kmer_length, split_kmers);
//                split_kmers = rolling_kmer_extract(current_contig, kmer_length, split_kmers);
//                split_kmers = rolling_kmer_extract(current_contig, kmer_length, split_kmers);
                split_kmers = rolling_kmer_bitvector(current_contig, kmer_length, split_kmers);
//                std::cout << "This dict has size: " << split_kmers.size() << std::endl;
                if (split_kmers.size() == 0) {
                #pragma omp critical
                    {
                        errors.push_back(std::runtime_error("No sequence found in " + fasta_path[sample_idx]));
                        interrupt = true;
                    }
                }
            }

            // creating output file
            std::ofstream os(names[sample_idx] + ".skf", std::ios::binary);
            // creating archive

            cereal::BinaryOutputArchive oarchive(os);

            MyBase cereal_split_kmers;
            cereal_split_kmers.cereal_dict = split_kmers;
            oarchive(cereal_split_kmers);
//            {
//                cereal::BinaryInputArchive iarchive(os);
//                cereal_class getting_back;
//                getting_back.cereal_dict;
//                iarchive(getting_back);
//            }

//                out_dict = getting_back;
//            }
//            for (const auto& iterator: out_dict) {
//                std::cout << iterator.first << ": " << iterator.second << std::endl;
//            }

            if (!interrupt) {
                // TODO: use std::move here rather than copy?
                kmer_dicts[sample_idx] = split_kmers;
                kseq_destroy(seq);
                gzclose(fp);
            }
        }
    }
    // Check for errors
    if (interrupt) {
        for (auto i = errors.cbegin(); i != errors.cend(); ++i) {
            std::cout << i->what() << std::endl;
        }
        if (errors.size()) {
            throw std::runtime_error("Errors while creating k-mer dictionaries!");
        } else {
            throw pybind11::error_already_set();
        }
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
//    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    std::cout << elapsed_seconds.count() << std::endl;
//    std::cout << "Got k-mers! " << std::endl;
    return kmer_dicts;
}

//TODO: map instead of unordered map
std::unordered_map<uint64_t, uint64_t> change_type(robin_hood::unordered_map<uint64_t, uint64_t> rh_dict) {
    std::unordered_map<uint64_t, uint64_t> basic_map;
    for (const auto& x: rh_dict) {
         basic_map[x.first] = x.second;
    }
    return basic_map;
}

vec_dict_bits ska_fasta(const std::vector<std::string>& isolate_vector, const std::vector<std::string>& isolate_names, int kmerLength) {
//std::vector<std::unordered_map<uint64_t, uint64_t>> ska_fasta(const std::vector<std::string>& isolate_vector, const std::vector<std::string>& isolate_names, int kmerLength) {

// kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
    vec_dict_bits kmer_dicts = get_kmers(isolate_vector, isolate_names, kmerLength);
//    std::vector<std::unordered_map<uint64_t, uint64_t>> kmer_dicts = get_kmers(isolate_vector, isolate_names, kmerLength);
    // TODO: pass kmer_dicts as reference to not copy everything over all the time
    return kmer_dicts;
}



int run_ska(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
////   from sketchlib: https://github.com/bacpop/pp-sketchlib/blob/master/src/api.cpp

    vec_dict_bits kmer_dicts = ska_fasta(isolate_paths, isolate_names, kmerLength);
    return 0;
}
