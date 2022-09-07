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
#include "progressbar.hpp"
#include <bitset>
#include <tuple>

//using namespace boost::filesystem;
KSEQ_INIT(gzFile, gzread)

// from Simon Hariss's original code
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
    std::string blub;
    uint64_t res = ~mer;
    std::cout << "Res: " << res << " , Mer: " << mer << std::endl;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    blub = (res >> (2 * (32 - kmerSize)));
    std::cout << (res >> (2 * (32 - kmerSize))) << ", " << blub << std::endl;
    return (res >> (2 * (32 - kmerSize)));
}

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
    for (int i1 = 0; i1 <= sequence.length() - k; i1++) {
//          TODO: rolling k-mers
        current_kmer = sequence.substr(i1, n) + sequence.substr(i4, n);
        split_kmer_base = sequence[i3];
        dict[current_kmer] = split_kmer_base;
        std::cout << current_kmer << ": " << split_kmer_base << std::endl;
        i3++;
        i4++;
    }
    return dict;
}

robin_hood::unordered_map<std::string, char> rolling_kmer_extract(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict) {
    std::string temp;
    std::string current_kmer;
    std::string split_kmer;
    char base;
    for (int i = k-1; i < sequence.length();  i++) {
        if (temp != "") {
            current_kmer = temp.substr(1, k-1) + sequence[i];
        }
        else {
            current_kmer = sequence.substr(0, k);
        }
//        base = extract_middlebase(current_kmer);
        base = current_kmer[(k/2)];
        temp = current_kmer;
        split_kmer = current_kmer.erase((k/2), 1);
        dict[split_kmer] = base;
        std::cout << current_kmer << ": " << base << std::endl;
    }
    return dict;
}

robin_hood::unordered_map<std::string, char> rolling_kmer_iterators(std::string sequence, int k, robin_hood::unordered_map<std::string, char> dict) {
    std::cout << sequence << std::endl;
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
vec_dict get_kmers(const std::vector< std::string>& fasta_path, const std::vector< std::string>& names, int kmer_length)
{
    progressbar bar1(fasta_path.size());
    std::cout << "Creating k-mer dictionary: " << std::endl;
    kmer_length = check_kmer_length(kmer_length);
    vec_dict kmer_dicts;
    kmer_dicts.resize(fasta_path.size());
// read in fasta files with kseq.h
    auto start = std::chrono::steady_clock::now();
    bool interrupt = false;
    std::vector<std::runtime_error> errors;


//#pragma omp parallel for // num_threads(2) // #pragma omp parallel is for paralyzing this loops
    for (int sample_idx = 0; sample_idx < fasta_path.size(); ++sample_idx) {
        #pragma omp critical
        bar1.update();
        if (interrupt || PyErr_CheckSignals() != 0) {
            interrupt = true;
        } else {
            //#pragma omp critical
            robin_hood::unordered_map<std::string, char> split_kmers;
            gzFile fp;
            kseq_t *seq;
            int number_of_seqs = 0;
//          TODO: change filesystem to boost
//          boost::filesystem::exists(fasta_path[sample_idx]);
            std::__fs::filesystem::path p(fasta_path[sample_idx]);
            if (!std::__fs::filesystem::exists(p)) {
                throw std::runtime_error("The given file does not exist!");
            }
            fp = gzopen(fasta_path[sample_idx].c_str(), "r");
            seq = kseq_init(fp);
            //open current kmer file to write to string
            while ((kseq_read(seq)) >= 0) {
                std::string current_contig = seq->seq.s;
                ++number_of_seqs;
//
//                split_kmers = rolling_kmer_iterators(current_contig, kmer_length, split_kmers);
//                split_kmers = rolling_kmer_extract(current_contig, kmer_length, split_kmers);
//                split_kmers = old_kmer_approach(current_contig, kmer_length, split_kmers);

                if (split_kmers.size() == 0) {
                #pragma omp critical
                    {
                        errors.push_back(std::runtime_error("No sequence found in " + fasta_path[sample_idx]));
                        interrupt = true;
                    }
                }
            }
            // print k-mer dict to file:
            std::ofstream current_kmer_file;
            current_kmer_file.open(names[sample_idx] + ".skf");
            std::stringstream bitstringstream;
            std::string bitstring = bitstringstream.str();

            for (const auto& x: split_kmers) {
                bitstringstream << x.first << ":" << x.second;
            }
            ascii_bitstring(bitstring);
            current_kmer_file << bitstring;
            current_kmer_file.close();

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
//    std::cout << "This is done! " << kmer_dicts.size() << std::endl;
    return kmer_dicts;
}


vec_dict ska_fasta(const std::vector<std::string>& isolate_vector, const std::vector<std::string>& isolate_names, int kmerLength) {
    // kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
    vec_dict kmer_dicts = get_kmers(isolate_vector, isolate_names, kmerLength);
    // TODO: pass kmer_dicts as reference to not copy everything over all the time
    return kmer_dicts;
}

uint64_t to_binary(std::string& current_kmer, int& length) {
    // convert k-mer to bitvector
    uint64_t packed_int = 0;
    for (auto it = current_kmer.cbegin(); it != current_kmer.cend(); ++it) {
        packed_int = packed_int << 2;
        switch (*(it)) {
            case 'A':
                packed_int += 0; // optimise out
                break;
            case 'C':
                packed_int += 1;
                break;
            case 'G':
                packed_int += 2;
                break;
            case 'T':
                packed_int += 3;
                break;
        }
    }
    std::bitset<64> x(packed_int);
    std::cout << x << " and " << packed_int << std::endl;
    return packed_int;
}



int run_ska(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
////   from sketchlib: https://github.com/bacpop/pp-sketchlib/blob/master/src/api.cpp
//
    std::string s = "ACTGAATC";
    std::cout << s << std::endl;

    robin_hood::unordered_map<std::string, char> kmer_dicts;
    std::vector<char> bases = {'A', 'C', 'G', 'T', '-'};
    create_bitvector_value(bases);

    ReverseComp64(01320031, 8);
//    vec_dict kmer_dicts = ska_fasta(isolate_paths, isolate_names, kmerLength);
//    robin_hood::unordered_map<std::string, std::vector<char>> all_dict;
//    all_dict = create_one_large_dictionary(kmer_dicts, kmerLength);
////    std::cout << "size of all_dict: " << all_dict.size() << std::endl;
//    std::ofstream outFile;
//    outFile.open("variant_alignment.aln");
//// printing create_one_large_dictionary in parallel
////    #pragma omp parallel for
//    for (int isolate_num = 0; isolate_num < isolate_paths.size(); isolate_num++)
//    {
//        outFile << ">" << isolate_names[isolate_num] << std::endl;
//        for (const auto& kmer: all_dict)
//        {
//            outFile << kmer.second[isolate_num];
//        }
//        outFile << "\n";
//    }
//    outFile.close();







//    for (int i1 = 0; i1 <= s.length() - kmerLength; i1++) {
////          TODO: rolling k-mers
//
//        int i3 = (kmerLength - 1) / 2;
////              int i2 = i3 - 1;
//        int i4 = i3 + 1;
////              int i5 = kmer_length - 1;
//        int n = (kmerLength - 1) / 2;
//        char split_kmer_base;
//        std::string current_kmer = s.substr(i1, n) + s.substr(i4, n);
//        std::string temp = "";
//    }






//        current_kmer = s.substr(i1, n) + s.substr(i4, n);
//        split_kmer_base = s[i3];
////        split_kmers[current_kmer] = split_kmer_base;
//        i3++;
//        i4++;
//        std::cout << current_kmer << std::endl;
//    }
    return 0;
}
