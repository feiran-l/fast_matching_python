#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "calc_geodesics.h"
#include <pybind11/embed.h>
#include <signal.h>

namespace py = pybind11;



void signalHandler(int signum) {
    exit(signum);
}


vector<vector<double> > fast_marching(string data_dir, bool verbose) {
    vector<vector<double> > res = calc_geodesic_matrix(data_dir, verbose);
    return res;
}


PYBIND11_MODULE(fast_marching, m) {
    py::bind_vector<vector<vector<double> > > (m, "FloatVector2D");
    m.doc() = "fast marching for mesh geodesics"; 

    signal(SIGINT, signalHandler);

    m.def("fast_marching", &fast_marching, "calculate pairwise geodesics of a mesh", 
          py::arg("data_dir"), py::arg("verbose"));
}
