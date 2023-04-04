//
// Created by Johanna Helene von Wachsmann on 28/07/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include "robin_hood.h"
//#include <cereal/archives/binary.hpp>
#include "progressbar.hpp"
#include "cereal/include/cereal/archives/binary.hpp"
#include "cereal/include/cereal/types/vector.hpp"
#include "cereal/include/cereal/types/string.hpp"
#include "eigen/Eigen/Eigen"
#include "ska.hpp"

static const char look_up_table2[16] = {'A','C', 'G', 'T', 'M', 'R', 'W', 'S', //
                                         'Y', 'K', 'V', 'H', 'D', 'B', 'N', 'X'// }

};

//
std::string print_kmers2(int kmer_length, uint64_t kmer)
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

std::vector<int64_t>inarchiving(std::string file, std::vector<int64_t> dst)
{
//    std::cout << "in archinving" << std::endl;
    std::ifstream os(file, std::ios::binary);
    cereal::BinaryInputArchive iarchive(os); // Create an input archive
    iarchive(dst);
    return dst;
}

robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> ska_align(std::vector<robin_hood::unordered_map<uint64_t, uint8_t>>& kmer_dicts, int kmerLength)
{
    robin_hood::unordered_map<int, int> base_count;
    robin_hood::unordered_map<std::uint64_t, std::vector<uint8_t>> alignemnt_dict;
    robin_hood::unordered_map<std::uint64_t, uint8_t> merge_dict;
    double percentage = 0.0;

    // merge all the keys together into a dictionary
    for (auto & kmer_dict : kmer_dicts)
    {
        merge_dict.insert(kmer_dict.begin(), kmer_dict.end());
    }
    // the counter is used to count in how many isolates the current k-mer is present. The k-mer is discarded if a
    // certain threshold isn't reached

    // iterate through the merged dictionary's keys (aka all k-mers) to check if they are also present in the
    // individual dictionaries

    for (auto& key: merge_dict)
    {
        double abundance = 0.0;
        bool variation = false;
        char prev = 'J';
        std::vector<uint8_t> bases;
        for (auto& kmer_dict : kmer_dicts)
        {
            if (kmer_dict.count(key.first) != 0)
            {
                abundance++;
                bases.push_back(kmer_dict[key.first]);
                if (prev != 'J' && prev != kmer_dict[key.first])
                {
                    variation = true;
                }
                prev = kmer_dict[key.first];
            }
            else
            {
                // if not present in this isolate we add a placeholder
                bases.push_back('-');
                variation = true;
                prev = '-';
            }
        }
        percentage = abundance/kmer_dicts.size();
    // TODO: change percentage

// filtering if kmer in too little isolates it is discarded, additionally if all the bases are the same because we only
// want the variant and not the same things.
        if (percentage >= .9 && variation)
        {
            alignemnt_dict[key.first] = bases;
        }
    bases.clear();
}
return alignemnt_dict;
}


void run_ska_align(const std::vector<std::string>& skf_paths, const std::vector<std::string>& isolate_names, int kmerLength, std::string output_directory, std::string cluster_name)
{
//    TODO: parallel
    std::vector<robin_hood::unordered_map<uint64_t, uint8_t>> all_kmers;
    all_kmers.resize(skf_paths.size());
    int counter = 0;
    for (std::string file : skf_paths)
    {
//        std::cout << "Read in path: " << file << std::endl;
        robin_hood::unordered_map<uint64_t, uint8_t> current_isolate_dict;
        // cereal needs to end in return to work
        std::vector<int64_t> dst;
        dst = inarchiving(file, dst);


//        std::vector<int> dst;
//        std::ifstream instream(file, std::ios::binary);
//        cereal::BinaryInputArchive iarchive(instream); // Create an input archive
//        iarchive(dst);

        for (int i = 0; i < dst.size(); i += 2)
        {
            uint8_t value = look_up_table2[dst[i+1]];
            current_isolate_dict[dst[i]] = value;
        }
//        std::cout << "dict: ";
//        for (auto elem : current_isolate_dict)
//        {
//            std::cout << print_kmers2(kmerLength, elem.first) << " " << elem.second << ", ";
////            std::cout << elem.first << " " << elem.second << ", ";
//        }
//        all_kmers.push_back(current_isolate_dict);
//        std::cout << "counter: " << counter << std::endl;
        all_kmers[counter] = current_isolate_dict;

        counter++;
    }
    robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> variant_alignemnt;
    variant_alignemnt = ska_align(all_kmers, kmerLength);

//    print variant alignment value
//    std::cout << "Variant alignemnt Size: " << isolate_names.size() << ", " << variant_alignemnt.size() << std::endl;
//    for (int i = 0; i < isolate_names.size(); i++) {
//        for (auto &it: variant_alignemnt) {
////            std::cout << "Variant alignemnt: " << i << " " << it.second[i] << std::endl;
//        }
//    }

    std::cout << "Done with ska_align!" << std::endl;
    //  TODO: make vectors for each file
    //  TODO: check if isolate names and vector order are correct
    //  TODO: cerealize the output to file as chars

    std::ofstream myfile;
    std::string path = output_directory + cluster_name + ".aln";
    myfile.open (path);
//    std::cout << "Printing alignments to file to: " << path << std::endl;
    for (int i = 0; i < isolate_names.size(); i++)
    {
        myfile << ">" + isolate_names[i] + "\n";
        for (auto& it: variant_alignemnt)
        {
            myfile << it.second[i];
        }
        myfile << "\n";
    }
    myfile.close();
}

std::vector<std::string>compare_kmers(std::string file, int k)
{
//    std::cout << "Read in path: " << file << std::endl;
    std::vector<std::string> return_vector;
    std::vector<int64_t> dst;

    dst = inarchiving(file, dst);


    for (int i = 0; i < dst.size()-1; i += 2)
    {
        char blub = look_up_table2[dst[i+1]];
        std::string temp;
        temp = blub;
        std::string s = print_kmers2(k, dst[i]) + ":" + temp;
        return_vector.push_back(s);
//        std::cout << "string: " << s << std::endl;
    }
    return return_vector;
}

