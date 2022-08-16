//
// Created by Johanna Helene von Wachsmann on 14/06/2022.
//
#include <fstream>
#include <vector>
#include <iostream>
#include <zlib.h>
#include <chrono>
#include "robin_hood.h"
#include "ska.hpp"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

int check_kmer_length(int length) {
    bool even = length % 2 == 0;
    if (even) {
        std::cout << "K-mers have to be odd in length so k-mer length was changed to " << length - 1 << "!"
                  << std::endl;
    }
    return length;
}

//TODO: robin_hood::unordered_node_map better than robin_hood::unordered_map?
vec_dict get_kmers(const std::vector< std::string>& fasta_path, int kmer_length)
{
    kmer_length = check_kmer_length(kmer_length);
    char split_kmer_base;
    std::string current_kmer;
    vec_dict kmer_dicts;
    kmer_dicts.reserve(fasta_path.size());

// read in fasta files with kseq.h
    auto start = std::chrono::steady_clock::now();
    #pragma omp parallel for // #pragma omp parallel is for paralyzing this loops
    for (auto name_it = fasta_path.begin(); name_it != fasta_path.end(); name_it++) {
        std::cout << "file path " << " " << name_it->c_str() << std::endl;
        robin_hood::unordered_map<std::string, char> split_kmers;
        int counter1 = 0;
        gzFile fp;
        kseq_t *seq;
        int number_of_seqs = 0;
        fp = gzopen(name_it->c_str(), "r");
        seq = kseq_init(fp);
        while ((kseq_read(seq)) >= 0){
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
            for (int i1 = 0; i1 <= current_contig.length() - kmer_length; i1++)
            {
                // create kmer for dictionary
//          TODO: rolling k-mers
//                std::cout << "current contig " << current_contig << " " << i1 << std::endl;
                current_kmer = current_contig.substr(i1, n) + current_contig.substr(i4, n);
                split_kmer_base = current_contig[i3];
                split_kmers[current_kmer] = split_kmer_base;
//                std::cout << "current k-mer " << current_kmer << ":" << split_kmer_base << std::endl;
                i3++; i4++; //i1++;
            }
        }
        //    print dictionary
//        for (const auto& x: split_kmers)
//        {
//            std::cout << "dictionary size: " << split_kmers.size() << " " << x.first << " " << x.second << std::endl;
            // printf("isolate dict no: %d\tdict size: %lu isolate: %s %s\n", i, kmer_dicts[i].size(), x.first.c_str(), x.second);
//        }
        kmer_dicts.push_back(split_kmers);
//        printf("Total n. of sequences: %d\n", number_of_seqs);
        kseq_destroy(seq);
        gzclose(fp);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
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

vec_dict ska_fasta(const std::vector<std::string>& isolate_vector, int kmerLength) {
    // kmer_dicts is a vector of length n containing n dicts; each dict holds all the split k-mers of the current isolate
    vec_dict kmer_dicts = get_kmers(isolate_vector, kmerLength);
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

int run_ska(const std::vector<std::string>& isolate_paths, const std::vector<std::string>& isolate_names, int kmerLength) {
    vec_dict kmer_dicts = ska_fasta(isolate_paths, kmerLength);
    robin_hood::unordered_map<std::string, std::vector<char>> all_dict;
    all_dict = create_one_large_dictionary(kmer_dicts, kmerLength);
    std::cout << "size of all_dict: " << all_dict.size() << std::endl;
    std::ofstream outFile;
    outFile.open("/Users/wachsmannj/Documents/SKA2/tests/variant_alignment.aln");
//// printing create_one_large_dictionary in parallel
    #pragma omp parallel for
    for (int isolate_num = 0; isolate_num < isolate_paths.size(); isolate_num++)
    {
        outFile << ">" << isolate_names[isolate_num] << std::endl;
        for (const auto& kmer: all_dict)
        {
            outFile << kmer.second[isolate_num];
        }
        outFile << "\n";
    }
    outFile.close();
    std::cout << "Done with alignment file!  " << std::endl;
    return 0;
}