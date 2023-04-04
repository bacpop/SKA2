/*
 * ska_bindings.cpp
 * Python bindings for ska
 *
 */

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ska.hpp"
namespace py = pybind11;

PYBIND11_MODULE(ska_cpp, m) {
  m.doc() = "SKA functions";

  // Exported functions
  m.def("run_ska_fasta", &run_ska_fasta, "Runs ska fasta",
        py::arg("path_strings"),
        py::arg("name_strings"),
        py::arg("kmer_length"),
        py::arg("output"));
  m.def("to_binary", &to_binary, "Turns string into binary representation",
        py::arg("string"),
        py::arg("kmer_length"));
  m.def("ReverseComp64", &ReverseComp64, "calculates canonical kmer",
        py::arg("kmer"),
        py::arg("kmer_length"));
  m.def("check_for_N", &check_for_N, "checks for Ns in sequence",
        py::arg("split"),
        py::arg("position"));
//  m.def("rolling_kmer_bitvector", &rolling_kmer_bitvector, "created integer k-mers",
//        py::arg("sequence"),
//        py::arg("kmer"),
//        py::arg("dictionary"));
  m.def("change_type", &change_type, "change type of dict for testing",
        py::arg("dictionary"));
  m.def("change_type2", &change_type2, "change type of dict for testing",
        py::arg("dictionary"));
  m.def("testing_run_ska", &testing_run_ska, "test ska fasta",
        py::arg("s"),
        py::arg("k"),
        py::arg("kmer_or_base"));
  m.def("run_ska_align", &run_ska_align, "Runs ska align",
        py::arg("skf_paths"),
        py::arg("isolate_names"),
        py::arg("kmerLength"),
        py::arg("output"),
        py::arg("cluster_name"));
  m.def("get_kmers", &get_kmers, "creates k-mers",
        py::arg("fasta_path"),
        py::arg("names"),
        py::arg("kmer_length"),
        py::arg("output"));
//  m.def("testing_run_ska_without_reverse_complement", &testing_run_ska_without_reverse_complement, "test_kmers",
//        py::arg("s"),
//        py::arg("k"),
//        py::arg("path"),
//        py::arg("file_name"));
//  m.def("testing_rolling_kmer_bitvector", &testing_rolling_kmer_bitvector, "created integer k-mers",
//        py::arg("sequence"),
//        py::arg("kmer"),
//        py::arg("dict"));
//  m.def("testing_run_ska_align_without_reverse_complement", &testing_run_ska_align_without_reverse_complement, "test_align",
//        py::arg("skf_files"),
//        py::arg("names"),
//        py::arg("k"));
  m.def("test_ambiguity_bases", &test_ambiguity_bases, "test_ambiguity_bases",
          py::arg("old_base"),
          py::arg("split_base"));
  m.def("test_check_for_Ns", &test_check_for_Ns, "test_check_for_Ns",
        py::arg("sequence"),
        py::arg("k"));
  m.def("new_rolling_kmer_bitvector", &new_rolling_kmer_bitvector, "new_rolling_kmer_bitvector",
        py::arg("sequence"),
        py::arg("k"),
        py::arg("dict"));
  m.def("test_new_rolling_kmer_binary", &test_new_rolling_kmer_binary, "test_new_rolling_kmer_binary",
        py::arg("sequence"),
        py::arg("k"));
  m.def("compare_kmers", &compare_kmers, "compare_kmers",
        py::arg("path"),
        py::arg("k"));
  m.def("print_kmers2", &print_kmers2, "print_kmers2",
        py::arg("k"),
        py::arg("kmer"));
  //TODO: https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-properly-handle-ctrl-c-in-long-running-functions



//  m.def("run_ska_align", &run_ska_align, "Runs ska fasta and align",
//        py::arg("fasta_strings"),
//        py::arg("kmer_length") = 7);
  // Example args
  /*
        py::arg("db_name"), py::arg("samples"), py::arg("files"),
        py::arg("klist"), py::arg("sketch_size"),
        py::arg("codon_phased") = false, py::arg("calc_random") = true,
        py::arg("use_rc") = true, py::arg("min_count") = 0,
        py::arg("exact") = false, py::arg("num_threads") = 1,
        py::arg("use_gpu") = false, py::arg("device_id") = 0);
  */

  m.attr("version") = VERSION_INFO;
}

