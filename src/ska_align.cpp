//
// Created by Johanna Helene von Wachsmann on 28/07/2022.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include "robin_hood.h"
#include <cereal/archives/binary.hpp>
#include "progressbar.hpp"
#include "cereal/include/cereal/archives/binary.hpp"
#include "cereal/include/cereal/types/vector.hpp"
#include "cereal/include/cereal/types/string.hpp"
#include "eigen/Eigen/Eigen"
#include "ska.hpp"
#include <eigen/Eigen/StdVector>
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Core>

static const char look_up_table2[16] = {'A','C', 'G', 'T', 'M', 'R', 'W', 'S', //
                                         'Y', 'K', 'V', 'H', 'D', 'B', 'N', 'X'// }

};

robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> ska_align(std::vector<robin_hood::unordered_map<uint64_t, uint8_t>>& kmer_dicts, int kmerLength)
{
    robin_hood::unordered_map<int, int> base_count;
    robin_hood::unordered_map<std::uint64_t, std::vector<uint8_t>> alignemnt_dict;
    robin_hood::unordered_map<std::uint64_t, uint8_t> merge_dict;
    double percentage = 0.0;

    // merge all the keys together into a dictionary
    for (auto & kmer_dict : kmer_dicts) {
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
                if (prev != 'J' && prev != kmer_dict[key.first]) {
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
    if (percentage > .9 && variation)
    {
        alignemnt_dict[key.first] = bases;
    }
    bases.clear();
}
return alignemnt_dict;
}


void run_ska_align(const std::vector<std::string>& skf_paths, const std::vector<std::string>& isolate_names, int kmerLength, std::string output_directory, std::string cluster_name) {
//    TODO: parallel
    std::vector<robin_hood::unordered_map<uint64_t, uint8_t>> all_kmers;
    for (std::string file : skf_paths) {
        robin_hood::unordered_map<uint64_t, uint8_t> current_isolate_dict;
        std::vector<int> dst;
        std::ifstream instream(file, std::ios::binary);
        cereal::BinaryInputArchive iarchive(instream); // Create an input archive
        iarchive(dst);

        for (int i = 0; i < dst.size(); i += 2) {
            uint8_t value = look_up_table2[dst[i+1]];
            current_isolate_dict[dst[i]] = value;
        }
        all_kmers.push_back(current_isolate_dict);
    }
    robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> variant_alignemnt;
    variant_alignemnt = ska_align(all_kmers, kmerLength);

    std::cout << "Done with ska_align!" << std::endl;
    //  TODO: make vectors for each file
    //  TODO: check if isolate names and vector order are correct
    //  TODO: cerealize the output to file as chars

    std::ofstream myfile;
    std::string path = output_directory + cluster_name + ".aln";
    myfile.open (path);
    std::cout << "Printing alignments to file to: " << path << std::endl;
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

void testing_run_ska_align_without_reverse_complement(std::vector<std::string> skf_paths, std::vector<std::string> isolate_names)
{
    std::vector<robin_hood::unordered_map<uint64_t, uint8_t>> all_kmers;
    for (std::string file : skf_paths)
    {
        robin_hood::unordered_map<uint64_t, uint8_t> current_isolate_dict;
        std::vector<int> dst;
        std::ifstream instream(file, std::ios::binary);
        cereal::BinaryInputArchive iarchive(instream); // Create an input archive
        iarchive(dst);
        for (int i = 0; i < dst.size(); i += 2)
        {
            uint8_t value = look_up_table2[dst[i+1]];
            current_isolate_dict[dst[i]] = value;
        }
        all_kmers.push_back(current_isolate_dict);
    }
    robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> variant_alignment;
    variant_alignment = ska_align(all_kmers, 31);
    std::cout << "Done with ska_align!" << std::endl;
    //  TODO: make vectors for each file
    //  TODO: check if isolate names and vector order are correct
    //  TODO: cerealize the output to file as chars
    std::ofstream myfile;
    std::string path = "/Users/wachsmannj/Documents/test_SKA2/test_all_without_reverse_complement/variant_alignment.aln";
    myfile.open (path);

    for (int i = 0; i < isolate_names.size(); i++)
    {
        std::cout << "printing alignments to file: " << isolate_names[i] << "!" << std::endl;
        myfile << ">" + isolate_names[i] + "\n";
        for (auto& it: variant_alignment)
        {
            myfile << it.second[i];
        }
        myfile << "\n";
    }
    myfile.close();
}
