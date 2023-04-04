////
//// Created by johanna on 15/03/23.
////
//
//void testing_run_ska_without_reverse_complement(std::vector<std::string> seqs, int k, std::string path, std::string file_name)
//{
////    normal run ska just without considering the reverse complement
//    std::unordered_map<uint64_t, uint8_t> dict;
//
//    for (auto s: seqs)
//    {
//
//        uint64_t kmer_mask = 0;
//        for (int i = 0; i < k / 2; i++)
//        {
//            kmer_mask = kmer_mask << 2;
//            kmer_mask += 3;
//        }
//        int pos = 0;
//        bool new_base_N = true;
//        uint64_t new_base;
//        uint8_t sub_kmer_length = k / 2;
//        std::string current_kmer = s.substr(pos, k);
//        std::vector<int> check;
//        uint64_t b1;
//        uint64_t b2;
//        uint64_t m;
//        uint64_t smallest_canonical;
//        uint64_t mask = 3;
//        while (new_base_N && pos < s.length())
//        {
//            check = check_for_N(current_kmer, pos);
//            // no mpre Ns
//            if (check[0] == 0) {
//                pos = check[1]; //9
//                new_base = look_up_table[s[k + pos]];
//                new_base_N = false;
//                // convert to binary
//                std::string split1 = current_kmer.substr(0, k / 2);
//                std::string split2 = current_kmer.substr(k / 2 + 1, k / 2);
//                b1 = to_binary(split1, k / 2);
//                b2 = to_binary(split2, k / 2);
//
//                uint64_t bint = (b1 << ((k / 2) + 1)) | b2;
//                m = look_up_table[s[sub_kmer_length]];
//                smallest_canonical = bint;
//                int64_t value = m;
//
//                // update value if kmer already existed
//                //key not present in dictionary
//                uint8_t mid_val = static_cast<uint8_t>(value);
//                if (dict.find(smallest_canonical) == dict.end())
//                {
//                    uint8_t u8_value1 = value;
//                    dict[smallest_canonical] = u8_value1;
//                }
//                else
//                {
//                    // calculate offset of position in vector
//                    int64_t value1 = ambiguous_bases[(m * 15) + dict[smallest_canonical]];
//                    uint8_t u8_value1 = value1;
//                    dict[smallest_canonical] = u8_value1;
////                    //                    TODO: change back to ambiguity
////                    uint8_t u8_value = value;
////                    dict[smallest_canonical] = u8_value;
//                }
//            }
//                // more Ns
//            else {
//                pos = check[1];
//                current_kmer = s.substr(pos, k);
//                new_base_N = true;
//            }
//        }
////TODO: does N still = 4 here??
//        new_base = look_up_table[s[k + pos]]; //should be either 00,01,10,11
//        for (int new_base_pos = pos; new_base_pos < s.size() - k; ++new_base_pos)
//        {
//            new_base = look_up_table[s[k + new_base_pos]];
////            std::string new_base_print = std::bitset<4>(new_base).to_string();
//            if (new_base == 14) {
//                new_base_N = true;
//                new_base_pos = new_base_pos + k + 1;
//                current_kmer = s.substr(new_base_pos, k);
//                while (new_base_N && new_base_pos < s.length())
//                {
//                    check = check_for_N(current_kmer, new_base_pos);
//                    // no mpre Ns
//                    if (check[0] == 0)
//                    {
//                        new_base_pos = check[1] - 1; //because of loop?
//                        new_base_N = false;
//                    }
//                        // more Ns
//                    else
//                    {
//                        new_base_pos = check[1];
//                        current_kmer = s.substr(new_base_pos, k);
//                        new_base_N = true;
//                    }
//                }
//            }
//            else
//            {
//                b1 = b1 << 2;
//                b1 += m;
//                b1 = b1 & kmer_mask;
//                m = b2 >> (((k / 2) * 2) - 2); //does this work??
//                new_base = look_up_table[s[k + new_base_pos]];
//                b2 = b2 << 2;
//                b2 += new_base;
//                b2 = b2 & kmer_mask;
//                uint64_t bitstring = (b1 << (k / 2 + 1)) | b2;
//                smallest_canonical = bitstring;
//                int64_t value = m;
//                if (smallest_canonical < bitstring)
//                {
//                    value = ~m;
//                    value = value & mask;
//                }
//                else
//                {
//                    smallest_canonical = bitstring;
//                }
//                // update value if kmer already existed
//                //key not present in dictionary
//                if (dict.find(smallest_canonical) == dict.end())
//                {
//                    uint8_t u8_value = value;
//                    dict[smallest_canonical] = u8_value;
//                }
//                else
//                {
//                    // calculate offset of position in vector
//                    int64_t value1 = ambiguous_bases[(value * 15) + dict[smallest_canonical]];
//                    uint8_t u8_value1 = value1;
//                    dict[smallest_canonical] = u8_value1;
//
//////                    TODO: change back to ambiguity
////                    uint8_t u8_value = value;
////                    dict[smallest_canonical] = u8_value;
//                }
//            }
//        }
//        std::vector<int> test_vector;
//
//        for (auto const &pair: dict)
//        {
//            test_vector.push_back(pair.first);
//            test_vector.push_back(pair.second);
////            std::cout << "{" << pair.first << ": " << look_up_table2[pair.second] << "}\n";
//        }
//
//        std::ofstream os("/home/johanna/test_SKA2/test_all_without_reverse_complement/" + file_name + ".skf",std::ios::binary);
//        cereal::BinaryOutputArchive oarchive(os);
//        oarchive(test_vector);
//    }
//}
////
////robin_hood::unordered_map<uint64_t, uint8_t> rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> dict)
////{
////    uint64_t kmer_mask = 0;
////    for (int i = 0; i < k/2; i++)
////    {
////        kmer_mask = kmer_mask << 2;
////        kmer_mask += 3;
////    }
////    int pos = 0;
////    bool new_base_N = true;
////    uint64_t new_base;
////    uint8_t sub_kmer_length = k/2;
////    std::string current_kmer = sequence.substr(pos, k);
////    std::vector<int> check;
////    uint64_t b1;
////    uint64_t b2;
////    uint64_t m;
////    uint64_t smallest_canonical;
////    uint64_t mask = 3;
////    while (new_base_N && pos < sequence.length())
////    {
////        check = check_for_N(current_kmer, pos);
////        // no mpre Ns
////        if (check[0] == 0)
////        {
////            pos = check[1]; //9
////            new_base = look_up_table[sequence[k + pos]] ;
////            new_base_N = false;
////            // convert to binary
////            std::string split1 = current_kmer.substr(0, k/2);
////            std::string split2 = current_kmer.substr(k/2 + 1, k/2);
////            b1 = to_binary(split1, k/2);
////            b2 = to_binary(split2, k/2);
////
////            uint64_t bint = (b1 << ((k/2)+1)) | b2;
////            m = look_up_table[sequence[sub_kmer_length]];
////            smallest_canonical = ReverseComp64(bint, k-1);
////            int64_t value = m;
////
////            if (smallest_canonical < bint)
////            {
////                value = ~m;
////                value = value & mask;
////            }
////            else
////            {
////                smallest_canonical = bint;
////            }
////            // update value if kmer already existed
////            //key not present in dictionary
//////            uint8_t mid_val = static_cast<uint8_t>(value);
////
////            if (dict.find(smallest_canonical) == dict.end())
////            {
////                uint8_t u8_value1 = value;
//////                print_kmers(k, smallest_canonical);
////                dict[smallest_canonical] = u8_value1;
////            }
////            else
////            {
////                // calculate offset of position in vector
////                int64_t value1 = ambiguous_bases[(m * 15) + dict[smallest_canonical]];
////                uint8_t u8_value1 = value1;
//////                print_kmers(k, smallest_canonical);
////                dict[smallest_canonical] = u8_value1;
////
////////                    TODO: change back to ambiguity
//////                uint8_t u8_value = value;
//////                dict[smallest_canonical] = u8_value;
////            }
////        }
////        // more Ns
////        else
////        {
////            pos = check[1];
////            current_kmer = sequence.substr(pos, k);
////            new_base_N = true;
////        }
////    }
//////TODO: does N still = 4 here??
////    new_base = look_up_table[sequence[k + pos]]; //should be either 00,01,10,11
////    for (int new_base_pos = pos; new_base_pos < sequence.size() - k; ++new_base_pos)
////    {
//////        std::cout << sequence[k + new_base_pos] << std::endl;
////        new_base = look_up_table[sequence[k + new_base_pos]];
//////        std::cout << new_base << std::endl;
//////        std::string new_base_print = std::bitset<4>(new_base).to_string();
////        if(new_base == 14)
////        {
////            new_base_N = true;
////            new_base_pos = new_base_pos + k + 1;
////            current_kmer = sequence.substr(new_base_pos, k);
////            while (new_base_N && new_base_pos < sequence.length())
////            {
////                check = check_for_N(current_kmer, new_base_pos);
////                // no mpre Ns
////                if (check[0] == 0)
////                {
////                    new_base_pos = check[1]-1; //because of loop?
////                    new_base_N = false;
////                }
////                // more Ns
////                else
////                {
////                    new_base_pos = check[1];
////                    current_kmer = sequence.substr(new_base_pos, k);
////                    new_base_N = true;
////                }
////            }
////        }
////        else
////        {
////            b1 = b1 << 2;
////            b1 += m;
////            b1 = b1 & kmer_mask;
////            m = b2 >> (((k/2)*2)-2); //TODO: does this work??
////            new_base = look_up_table[sequence[k + new_base_pos]];
////            b2 = b2 << 2;
////            b2 += new_base;
////            b2 = b2 & kmer_mask;
////            uint64_t bitstring = (b1 << (k/2+1)) | b2;
////            smallest_canonical = ReverseComp64(bitstring, k/2+1);
////            int64_t value = m;
////            if (smallest_canonical < bitstring)
////            {
////                value = ~m;
////                value = value & mask;
////            }
////            else
////            {
////                smallest_canonical = bitstring;
////            }
////            // update value if kmer already existed
////            //key not present in dictionary
////            if (dict.find(smallest_canonical) == dict.end())
////            {
////                uint8_t u8_value = value;
//////                print_kmers(k, smallest_canonical);
////                dict[smallest_canonical] = u8_value;
////            }
////            else
////            {
////                // calculate offset of position in vector
////                int64_t value1 = ambiguous_bases[(value * 15) + dict[smallest_canonical]];
////                uint8_t u8_value1 = value1;
////                dict[smallest_canonical] = u8_value1;
////////                    TODO: change back to ambiguity
//////                uint8_t u8_value = value;
//////                dict[smallest_canonical] = u8_value;
////
////            }
////        }
////    }
////    return dict;
////}
//
//std::unordered_map<uint64_t, uint8_t> testing_rolling_kmer_bitvector(std::string& sequence, int k, std::unordered_map<uint64_t, uint8_t> dict)
//{
//    uint64_t kmer_mask = 0;
//    for (int i = 0; i < k/2; i++)
//    {
//        kmer_mask = kmer_mask << 2;
//        kmer_mask += 3;
//    }
//    int pos = 0;
//    bool new_base_N = true;
//    uint64_t new_base;
//    uint8_t sub_kmer_length = k/2;
//    std::string current_kmer = sequence.substr(pos, k);
//    std::vector<int> check;
//    uint64_t b1;
//    uint64_t b2;
//    uint64_t m;
//    uint64_t smallest_canonical;
//    uint64_t mask = 3;
//    while (new_base_N && pos < sequence.length())
//    {
//        check = check_for_N(current_kmer, pos);
//        // no mpre Ns
//        if (check[0] == 0)
//        {
//            pos = check[1]; //9
//            new_base = look_up_table[sequence[k + pos]] ;
//            new_base_N = false;
//            // convert to binary
//            std::string split1 = current_kmer.substr(0, k/2);
//            std::string split2 = current_kmer.substr(k/2 + 1, k/2);
//            b1 = to_binary(split1, k/2);
//            b2 = to_binary(split2, k/2);
//
//            uint64_t bint = (b1 << ((k/2)+1)) | b2;
//            m = look_up_table[sequence[sub_kmer_length]];
//            smallest_canonical = ReverseComp64(bint, k-1);
//            int64_t value = m;
//
//            if (smallest_canonical < bint)
//            {
//                value = ~m;
//                value = value & mask;
//            }
//            else
//            {
//                smallest_canonical = bint;
//            }
//            // update value if kmer already existed
//            //key not present in dictionary
////            uint8_t mid_val = static_cast<uint8_t>(value);
//
//            if (dict.find(smallest_canonical) == dict.end())
//            {
//                uint8_t u8_value1 = value;
////                print_kmers(k, smallest_canonical);
//                dict[smallest_canonical] = u8_value1;
//            }
//            else
//            {
//                // calculate offset of position in vector
//                int64_t value1 = ambiguous_bases[(m * 15) + dict[smallest_canonical]];
//                uint8_t u8_value1 = value1;
////                print_kmers(k, smallest_canonical);
//                dict[smallest_canonical] = u8_value1;
//
//////                    TODO: change back to ambiguity
////                uint8_t u8_value = value;
////                dict[smallest_canonical] = u8_value;
//            }
//        }
//            // more Ns
//        else
//        {
//            pos = check[1];
//            current_kmer = sequence.substr(pos, k);
//            new_base_N = true;
//        }
//    }
////TODO: does N still = 4 here??
//    new_base = look_up_table[sequence[k + pos]]; //should be either 00,01,10,11
//    for (int new_base_pos = pos; new_base_pos < sequence.size() - k; ++new_base_pos)
//    {
////        std::cout << sequence[k + new_base_pos] << std::endl;
//        new_base = look_up_table[sequence[k + new_base_pos]];
////        std::cout << new_base << std::endl;
////        std::string new_base_print = std::bitset<4>(new_base).to_string();
//        if(new_base == 14)
//        {
//            new_base_N = true;
//            new_base_pos = new_base_pos + k + 1;
//            current_kmer = sequence.substr(new_base_pos, k);
//            while (new_base_N && new_base_pos < sequence.length())
//            {
//                //    TODO: are Ns skipped when they are the new base to be added to the k-mer?
//                check = check_for_N(current_kmer, new_base_pos);
//                // no mpre Ns
//                if (check[0] == 0)
//                {
//                    new_base_pos = check[1]-1; //because of loop?
//                    new_base_N = false;
//                }
//                    // more Ns
//                else
//                {
//                    new_base_pos = check[1];
//                    current_kmer = sequence.substr(new_base_pos, k);
//                    new_base_N = true;
//                }
//            }
//        }
//        else
//        {
//            b1 = b1 << 2;
//            b1 += m;
//            b1 = b1 & kmer_mask;
//            m = b2 >> (((k/2)*2)-2); //TODO: does this work??
//            new_base = look_up_table[sequence[k + new_base_pos]];
//            b2 = b2 << 2;
//            b2 += new_base;
//            b2 = b2 & kmer_mask;
//            uint64_t bitstring = (b1 << (k/2+1)) | b2;
//            smallest_canonical = ReverseComp64(bitstring, k/2+1);
//            int64_t value = m;
//            if (smallest_canonical < bitstring)
//            {
//                value = ~m;
//                value = value & mask;
//            }
//            else
//            {
//                smallest_canonical = bitstring;
//            }
//            // update value if kmer already existed
//            //key not present in dictionary
//            if (dict.find(smallest_canonical) == dict.end())
//            {
//                uint8_t u8_value = value;
////                print_kmers(k, smallest_canonical);
//                dict[smallest_canonical] = u8_value;
//            }
//            else
//            {
//                // calculate offset of position in vector
//                int64_t value1 = ambiguous_bases[(value * 15) + dict[smallest_canonical]];
//                uint8_t u8_value1 = value1;
//                dict[smallest_canonical] = u8_value1;
//////                    TODO: change back to ambiguity
////                uint8_t u8_value = value;
////                dict[smallest_canonical] = u8_value;
//
//            }
//        }
//    }
//    return dict;
//}
//
//void get_kmers()
//{
//
//////            normal printing:
////            std::ofstream current_kmer_file;
////            current_kmer_file.open("/Users/wachsmannj/Documents/test_SKA2/integer_approach/" + names[sample_idx] + ".skf");
//////            std::cout << names[sample_idx] << ".skf :" << split_kmers.size() << std::endl;
////            std::stringstream bitstringstream;
////            for (const auto& x: split_kmers) {
////                bitstringstream << x.first << " " << x.second << " ";
////            }
////            current_kmer_file << bitstringstream.str();
////            current_kmer_file.close();
////            auto end1 = std::chrono::steady_clock::now();
////            std::chrono::duration<double> elapsed_seconds_kmers = end1-start;
////            std::cout << "elapsed time: " << elapsed_seconds_kmers.count() << "s\n";
//
////            // creating output file with cereal dictionary
////            std::ofstream os("/Users/wachsmannj/Documents/test_SKA2/integer_approach/testing/" +names[sample_idx] + ".skf", std::ios::binary);
////            cereal::BinaryOutputArchive oarchive(os);
////            MyBase cereal_split_kmers;
////            cereal_split_kmers.cereal_dict = split_kmers;
////            oarchive(cereal_split_kmers);
////            auto end1 = std::chrono::steady_clock::now();
////            std::chrono::duration<double> elapsed_seconds1 = end1-start;
//
////            std::ofstream os1("");
////            std::ofstream mystringfile;
////            mystringfile.open (output_directory + "/" + names[sample_idx] + "_string.skf");
//}
//
//
//robin_hood::unordered_map<uint64_t, uint8_t> rolling_kmer_bitvector(std::string& sequence, int k, robin_hood::unordered_map<uint64_t, uint8_t> dict)
//{
//    std::cout << "Rolling k-mers: " << sequence << std::endl;
//    uint64_t kmer_mask = 0;
//    for (int i = 0; i < k/2; i++)
//    {
//        kmer_mask = kmer_mask << 2;
//        kmer_mask += 3;
//    }
//    int pos = 0;
//    bool new_base_N = true;
//    uint64_t new_base;
//    uint8_t sub_kmer_length = k/2;
//    std::string current_kmer = sequence.substr(pos, k);
//    std::cout << "current k-mer: " << current_kmer << std::endl;
//    std::vector<int> check;
//    uint64_t b1;
//    uint64_t b2;
//    uint64_t m;
//    uint64_t smallest_canonical;
//    uint64_t mask = 3;
////    ---------------------------------------------------------------------------------
//    while (new_base_N && pos < sequence.length())
//    {
//        //    TODO: are Ns skipped when they are the new base to be added to the k-mer? NOOOOOOOOOOOOOOOO
//        std::cout << "At top!" << std::endl;
//        check = check_for_N(current_kmer, pos);
//        // no more Ns
//        if (check[0] == 0)
//        {
//            pos = check[1]; //9
//            new_base = look_up_table[sequence[k + pos]];
//            std::cout << "new base: " << look_up_table2[new_base] << std::endl;
//            if (look_up_table2[new_base] == 'N' or look_up_table2[new_base] == 'n')
//            {
//                pos += k;
//                std::cout << "True!" << std::endl;
//                std::cout << "new base is N so we have to skip over it" << std::endl;
//                new_base_N = true;
//                current_kmer = sequence.substr(k + pos + 1, k);
//                std::cout << "current_kmer: " << current_kmer << std::endl;
//
//                new_base = look_up_table[sequence[k + pos]];
//                std::cout << "new base: " << look_up_table2[new_base] << std::endl;
//                break;
//            }
//            else
//            {
//                new_base_N = false;
//            }
//
//            // convert to binary
//            std::string split1 = current_kmer.substr(0, k/2);
//
//            std::string split2 = current_kmer.substr(k/2 + 1, k/2);
//            b1 = to_binary(split1, k/2);
//            b2 = to_binary(split2, k/2);
//
//            uint64_t bint = (b1 << ((k/2)+1)) | b2; //binary_split_kmer
//            m = look_up_table[sequence[sub_kmer_length]]; //split_base
//            smallest_canonical = ReverseComp64(bint, k-1);
//            int64_t value = m; //prev_split_base
//
//            if (smallest_canonical < bint)
//            {
//                value = ~m;
//                value = value & mask;
//            }
//            else
//            {
//                smallest_canonical = bint;
//            }
//            // update value if kmer already existed
//            //key not present in dictionary
////            uint8_t mid_val = static_cast<uint8_t>(value);
//
//            if (dict.find(smallest_canonical) == dict.end())
//            {
//                uint8_t u8_value1 = value;
////                print_kmers(k, smallest_canonical);
//                dict[smallest_canonical] = u8_value1;
//            }
//            else
//            {
//                // calculate offset of position in vector
//                int64_t value1 = ambiguous_bases[(m * 15) + dict[smallest_canonical]];
//                uint8_t u8_value1 = value1;
////                print_kmers(k, smallest_canonical);
//                dict[smallest_canonical] = u8_value1;
//
//////                    TODO: change back to ambiguity
////                uint8_t u8_value = value;
////                dict[smallest_canonical] = u8_value;
//            }
//        }
//            // more Ns
//        else
//        {
//            pos = check[1];
//            current_kmer = sequence.substr(pos, k);
//            new_base_N = true;
//        }
//    }
////TODO: does N still = 4 here??
//    new_base = look_up_table[sequence[k + pos]]; //should be either 00,01,10,11
//    for (int new_base_pos = pos; new_base_pos < sequence.size() - k; ++new_base_pos)
//    {
////        std::cout << sequence[k + new_base_pos] << std::endl;
//        new_base = look_up_table[sequence[k + new_base_pos]];
////        std::cout << new_base << std::endl;
////        std::string new_base_print = std::bitset<4>(new_base).to_string();
//        if(new_base == 14)
//        {
//            new_base_N = true;
//            new_base_pos = new_base_pos + k + 1;
//            current_kmer = sequence.substr(new_base_pos, k);
//            while (new_base_N && new_base_pos < sequence.length())
//            {
//                check = check_for_N(current_kmer, new_base_pos);
//                // no mpre Ns
//                if (check[0] == 0)
//                {
//                    new_base_pos = check[1]-1; //because of loop?
//                    new_base_N = false;
//                }
//                    // more Ns
//                else
//                {
//                    new_base_pos = check[1];
//                    current_kmer = sequence.substr(new_base_pos, k);
//                    new_base_N = true;
//                }
//            }
//        }
//        else
//        {
//            b1 = b1 << 2;
//            b1 += m;
//            b1 = b1 & kmer_mask;
//            m = b2 >> (((k/2)*2)-2); //TODO: does this work??
//            new_base = look_up_table[sequence[k + new_base_pos]];
//            b2 = b2 << 2;
//            b2 += new_base;
//            b2 = b2 & kmer_mask;
//            uint64_t bitstring = (b1 << (k/2+1)) | b2;
//            smallest_canonical = ReverseComp64(bitstring, k/2+1);
//            int64_t value = m;
//            if (smallest_canonical < bitstring)
//            {
//                value = ~m;
//                value = value & mask;
//            }
//            else
//            {
//                smallest_canonical = bitstring;
//            }
//            // update value if kmer already existed
//            //key not present in dictionary
//            if (dict.find(smallest_canonical) == dict.end())
//            {
//                uint8_t u8_value = value;
////                print_kmers(k, smallest_canonical);
//                dict[smallest_canonical] = u8_value;
//            }
//            else
//            {
//                // calculate offset of position in vector
//                int64_t value1 = ambiguous_bases[(value * 15) + dict[smallest_canonical]];
//                uint8_t u8_value1 = value1;
//                dict[smallest_canonical] = u8_value1;
//////                    TODO: change back to ambiguity
////                uint8_t u8_value = value;
////                dict[smallest_canonical] = u8_value;
//
//            }
//        }
//    }
//    for (const auto& x: dict)
//    {
//        std::cout << "Split k-mer: " << print_kmers(k, x.first) << " and base: " << look_up_table2[x.second] << std::endl;
//    }
//    return dict;
//}

//
//void testing_run_ska_align_without_reverse_complement(std::vector<std::string> skf_paths, std::vector<std::string> isolate_names, int k)
//{
//    std::vector<robin_hood::unordered_map<uint64_t, uint8_t>> all_kmers;
//    for (std::string file : skf_paths)
//    {
//        std::cout << "Here: " << file << std::endl;
//        robin_hood::unordered_map<uint64_t, uint8_t> current_isolate_dict;
//        std::vector<int> dst;
//        std::ifstream instream(file, std::ios::binary);
//        cereal::BinaryInputArchive iarchive(instream); // Create an input archive
//        iarchive(dst);
//        for (int i = 0; i < dst.size(); i += 2)
//        {
//            uint8_t value = look_up_table2[dst[i+1]];
//            current_isolate_dict[dst[i]] = value;
//        }
//        all_kmers.push_back(current_isolate_dict);
//        std::cout << "Here2: " << all_kmers.size() << std::endl;
//    }
//    robin_hood::unordered_map<std::uint64_t , std::vector<uint8_t>> variant_alignment;
//    variant_alignment = ska_align(all_kmers, k);
//    std::cout << "Done with ska_align!" << std::endl;
//    //  TODO: make vectors for each file
//    //  TODO: check if isolate names and vector order are correct
//    //  TODO: cerealize the output to file as chars
//    std::ofstream myfile;
//    std::string path = "/home/johanna/test_SKA2/integer_approach/MA_results/variant_alignment.aln";
//    myfile.open (path);
//
//    for (int i = 0; i < isolate_names.size(); i++)
//    {
//        std::cout << "printing alignments to file: " << isolate_names[i] << "!" << std::endl;
//        myfile << ">" + isolate_names[i] + "\n";
//        for (auto& it: variant_alignment)
//        {
//            myfile << it.second[i];
//        }
//        myfile << "\n";
//    }
//    myfile.close();
//}
