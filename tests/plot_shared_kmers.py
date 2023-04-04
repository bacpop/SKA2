import ska_cpp
import pandas as pd
import collections
import struct
import numpy as np

def main():
    k = 31
    # unfiltered
    # should be list with paths to all skf files
    # skf_paths = "/home/johanna/test_SKA2/integer_approach/test_speed_memory/rfile.skf"
    # print("Cluster 2 (unfiltered)")
    # with open(skf_paths) as f:
    #     skf_files = [line.rstrip() for line in f]
    #
    # print(skf_files[0])
    # kmers_in_intersection = set(ska_cpp.compare_kmers(skf_files[0], k))
    #
    # for i in range(1, len(skf_files)):
    #     compare_list = set(ska_cpp.compare_kmers(skf_files[i], k))
    #     kmers_in_intersection = kmers_in_intersection.intersection(compare_list)
    #     print(len(kmers_in_intersection))

    # filtered
    skf_paths = "/home/johanna/test_SKA2/integer_approach/skf_cluster2_unfiltered/rfile.skf"
    print("Cluster 2 (unfiltered)")
    with open(skf_paths) as f:
        skf_files = [line.rstrip() for line in f]

    print(skf_files[0])
    print(len(skf_files))
    vector_of_dicts = {}
    for i in range(0, 2):#len(skf_files)):
        print("i", i)
        current_dict = {}
        items = ska_cpp.compare_kmers(skf_files[i], k)
        for item in items:
            item = item.split(":")
            current_dict[item[0]] = item[1]
        vector_of_dicts.append(current_dict)

    merge_vector = []
    counter_1 = 1
    for dic in vector_of_dicts:
        print("# of dict:", counter_1, "/", len(skf_files))
        counter_1 += 1
        merge_vector.append(list(dic.keys()))
    flat_merge_vector = [item for blub in merge_vector for item in blub]
    myset = set(flat_merge_vector)
    unique_merge_vector = list(myset)
    unique_merge_vector.sort()
    print("Done Merging")


    alignment_dict = {}

    # TODO: filter out:
    #     1. abundance above 0.5%
    #     2. if vector all the same bases (bc then there is no variance)
    for kmer in unique_merge_vector:
        abundance = 0
        variation = False
        previous = 'J'
        current_base_vector = []
        for isolate_dict in unfiltered_list:
            # print('isolate_dict', isolate_dict)
            if kmer in isolate_dict:
                # print('kmer', kmer)
                abundance += 1
                current_base_vector.append(isolate_dict[kmer])
                if previous != 'J' and previous != isolate_dict[kmer]:
                    variation = True
                previous = isolate_dict[kmer]
            else:
                current_base_vector.append("-")
                variation = True
                previous = '-'
        # print(abundance, len(sequences))
        if abundance / len(sequences) > 0.9 and variation:
            alignment_dict[kmer] = current_base_vector
    vector_of_variants = []

    for i in range(0, len(paths)):
        sm = []
        for kmer in alignment_dict:
            sm.append(alignment_dict[kmer][i])
        sm.sort()
        vector_of_variants.append(sm)

    c_vector = []

    for i in range(0, len(paths)):
        sm = []
        for j in range(0, len(ska_seqs[0])):
            # print(ska_seqs[i][j])
            sm.append(ska_seqs[i][j])
        sm.sort()
        c_vector.append(sm)

    # print("merge List", len(merge_list))
    # # for isolate in unfiltered_list:
    # #     print("unfiltered", len(isolate))
    # merge_list = list(merge_list)
    # # for kmer in merge_list:
    # for kmer_index in range(0, len(merge_list)):
    #     count = 0
    #     # for current_isolate in unfiltered_list:
    #     for current_isolate_index in range(0, len(unfiltered_list)):
    #         if merge_list[kmer_index] in unfiltered_list[current_isolate_index]:
    #             count += 1
    #     if (count/len(unfiltered_list) >= 0.5):
    #     # print(kmer, count, count/len(unfiltered_list) < 0.5)
    #     # if (count/len(unfiltered_list) < 0.5 or count/len(unfiltered_list) == 1):
    #         print("here kmer", merge_list[kmer_index], current_isolate_index)
    #         # for current_isolate in unfiltered_list:
    #         #     if kmer in current_isolate:
    #         #         current_isolate.remove(kmer)
    #         #         print("removing kmer", kmer)
    #
    # print("Done Filtering")
    #
    # # for isolate in unfiltered_list:
    # #     print("filtered", len(isolate))

    print("# of dicts:", len(unfiltered_list))
    kmers_in_intersection = set(unfiltered_list[0])
    for i in range(1, len(unfiltered_list)):
        compare_list = set(unfiltered_list[i])
        kmers_in_intersection = kmers_in_intersection.intersection(compare_list)
        print(len(kmers_in_intersection))
        print(kmers_in_intersection)


if __name__ == "__main__":
    main()