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
#include "api.hpp"
#include <pybind11/pybind11.h>
#include <filesystem>
#include <iostream>
#include "progressbar.hpp"
#include <bitset>

//using namespace boost::filesystem;
KSEQ_INIT(gzFile, gzread)

// from Simon Hariss original code
void ascii_bitstring(std::string & mybits){
//    std::cout << "string: " << mybits << std::endl;
    int myremainder=::fmod(int(mybits.length()),6);
//    std::cout << "myremainder: " << myremainder << std::endl;
    if (myremainder>0){
        for (int i = 0; i<(6-myremainder); ++i){
            mybits.push_back('0');
//            std::cout << "mybits: " << mybits << std::endl;
        }
    }
    for (std::string::size_type i = 0; i<mybits.length(); i+=6){
        assert((i+5)<mybits.length());
        mybits[i/6] = ((int(mybits[i])-'0')+((int(mybits[i+1])-'0')*2)+((int(mybits[i+2])-'0')*4)+((int(mybits[i+3])-'0')*8)+((int(mybits[i+4])-'0')*16)+((int(mybits[i+5])-'0')*32))+33;
//        std::cout << "Mybits: " << mybits[i/6] << std::endl;
        }

    mybits.erase(mybits.length()/6);
}

int check_kmer_length(int length) {
    bool even = length % 2 == 0;
    bool thirtyone = length > 31;
    if (even) {
        std::cout << "K-mers have to be odd in length so k-mer length was changed to " << length - 1 << "!"
                  << std::endl;
    }
    else if (thirtyone) {
        std::cout << "K-mers cannot be longer than 31!"
                  << std::endl;
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
//        std::cout << sample_idx << std::endl;
#pragma omp critical
        bar1.update();
        if (interrupt || PyErr_CheckSignals() != 0) {
            interrupt = true;
        } else {
            //#pragma omp critical
            //std::cout << "file path " << " " << fasta_path[sample_idx] << std::endl;
            robin_hood::unordered_map<std::string, char> split_kmers;
            int counter1 = 0;
            gzFile fp;
            kseq_t *seq;
            int number_of_seqs = 0;
//            std::cout << "File path: " << fasta_path[sample_idx] << std::endl;
//            TODO: change filesystem to boost
//            boost::filesystem::exists(fasta_path[sample_idx]);
            std::__fs::filesystem::path p(fasta_path[sample_idx]);
            if (!std::__fs::filesystem::exists(p)) {
//                std::cout << fasta_path[sample_idx] << std::endl;
                throw std::runtime_error("The given file does not exist!");
            }
            fp = gzopen(fasta_path[sample_idx].c_str(), "r");
            seq = kseq_init(fp);
            //open current kmer file to write to string

//            std::cout << "names: " << names[sample_idx] << std::endl;

            while ((kseq_read(seq)) >= 0) {
//            std::cout << counter1 << std::endl;
                counter1++;
                std::string current_contig = seq->seq.s;
//            printf("seq blub: %s\n", seq->seq.s);
                ++number_of_seqs;
                // initializing the iterators which are needed to create the kmers
//            int i1 = 0;
                int i3 = (kmer_length - 1) / 2;
//            int i2 = i3 - 1;
                int i4 = i3 + 1;
//            int i5 = kmer_length - 1;
                int n = (kmer_length - 1) / 2;

                char split_kmer_base;
                std::string current_kmer;
                for (int i1 = 0; i1 <= current_contig.length() - kmer_length; i1++) {
//          TODO: rolling k-mers
//                std::cout << "current contig " << current_contig << " " << i1 << std::endl;
                    current_kmer = current_contig.substr(i1, n) + current_contig.substr(i4, n);
                    split_kmer_base = current_contig[i3];
                    split_kmers[current_kmer] = split_kmer_base;
                    i3++;
                    i4++; //i1++;
                }
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
//            current_kmer_file << "SKA __version__  " << std::endl;
//            current_kmer_file << names[sample_idx] << std::endl;
//            current_kmer_file << kmer_length << std::endl;
//            for (const auto& x: split_kmers) {
//                current_kmer_file << x.first << ": " << x.second << std::endl;
//            }
//            current_kmer_file.close();

            std::stringstream bitstringstream;
            std::string bitstring = bitstringstream.str();

            for (const auto& x: split_kmers) {
//                current_kmer_file << x.first << ": " << x.second << std::endl;
                bitstringstream << x.first;
            }
//            std::string bitstring = bitstringstream.str();
            ascii_bitstring(bitstring);
            current_kmer_file << bitstring;

//            for (const auto& x: split_kmers) {
//                    current_kmer_file << x.first;
//                current_kmer_file << std::endl;
//            }

            current_kmer_file.close();

            if (!interrupt) {
                // TODO: use std::move here rather than copy?
                kmer_dicts[sample_idx] = split_kmers;
//        printf("Total n. of sequences: %d\n", number_of_seqs);
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

//    std::cout << kmer_dicts.size() << std::endl;
////    print dictionary
//    for (int j = 0; j < kmer_dicts.size(); j++)
//    {
////        std::cout << "something size: " << kmer_dicts[j].size() << std::endl;
//        for (const auto& x: kmer_dicts[j])
//        {
//            std::cout << "isolate dict no: " << j << ", dict size:" << kmer_dicts[j].size() << " isolate " << x.first << " "
//                      << x.second << std::endl;
//            // printf("isolate dict no: %d\tdict size: %lu isolate: %s %s\n", i, kmer_dicts[i].size(), x.first.c_str(), x.second);
//        }
////        std::cout << std::endl;
//    }
    return kmer_dicts;
}

//std::vector<char> create_bitvector_value(const std::vector<char>& bases, int dim) {
//TODO: turn vector of char into bitvector to be more space efficient?
////    bases = [G, C, G] and bases = [A, -, A]
//}

vec_dict ska_fasta(const std::vector<std::string>& isolate_vector, const std::vector<std::string>& isolate_names, int kmerLength) {
    // kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
    vec_dict kmer_dicts = get_kmers(isolate_vector, isolate_names, kmerLength);
//    vec_dict get_kmers(const std::vector< std::string>& isolates, int kmer_length, const std::vector< int>& contigCount);
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

    return kmer_dicts;
}

void to_bitset(std::string& current_kmer, int& length) {
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
    std::bitset<8> x(packed_int);
    std::cout << x << " and " << packed_int << std::endl;
}


int run_ska(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
//   from sketchlib: https://github.com/bacpop/pp-sketchlib/blob/master/src/api.cpp

//    vec_dict kmer_dicts = ska_fasta(isolate_paths, isolate_names, kmerLength);
//    robin_hood::unordered_map<std::string, std::vector<char>> all_dict;
//    all_dict = create_one_large_dictionary(kmer_dicts, kmerLength);
////    std::cout << "size of all_dict: " << all_dict.size() << std::endl;
//    std::ofstream outFile;
//    outFile.open("variant_alignment.aln");
////// printing create_one_large_dictionary in parallel
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
    std::cout << "Done with alignment file!  " << std::endl;
    std::string current_kmer = "ATTA";
    to_bitset(current_kmer, kmerLength);

    return 0;
}