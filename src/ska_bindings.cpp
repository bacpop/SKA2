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
  m.def("run_ska", &run_ska, "Runs ska fasta and align",
        py::arg("fasta_strings"),
        py::arg("contig_count"),
        py::arg("kmer_length") = 7);
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

