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
        py::arg("kmer_length") = 7);
  m.def("to_binary", &to_binary, "Turns string into binary representation",
        py::arg("string"),
        py::arg("kmer_length"));
  m.def("ReverseComp64", &ReverseComp64, "calculates canonical kmer",
        py::arg("kmer"),
        py::arg("kmer_length"));
  m.def("check_for_N", &check_for_N, "checks for Ns in sequence",
        py::arg("split"),
        py::arg("position"));
  m.def("rolling_kmer_bitvector", &rolling_kmer_bitvector, "created integer k-mers",
        py::arg("sequence"),
        py::arg("kmer"),
        py::arg("dictionary"));
  m.def("change_type", &change_type, "change type of dict for testing",
        py::arg("dictionary"));
  m.def("change_type2", &change_type2, "change type of dict for testing",
        py::arg("dictionary"));
  m.def("testing_run_ska", &testing_run_ska, "test ska fasta",
        py::arg("sequence"),
        py::arg("kmer_length"));
  m.def("run_ska_align", &run_ska_align, "Runs ska align",
        py::arg("skf_paths"),
        py::arg("isolate_names"),
        py::arg("kmerLength"));
  m.def("get_kmers", &get_kmers, "creates k-mers",
        py::arg("fasta_path"),
        py::arg("names"),
        py::arg("kmer_length"));


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

