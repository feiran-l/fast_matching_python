#include "calc_geodesics.h"
#include "che_off.h"
#include "progress_bar.h"



vector<vector<double> > calc_geodesic_matrix(const string data_path, bool verbose) {

    che *mesh = new che_off(data_path);
    size_t n_vertices = mesh->n_vertices();
    vector<vector<double> > res(n_vertices);
    progressbar bar(n_vertices);

    // loop through each vertex
#pragma omp parallel for
    for (size_t i = 0; i < n_vertices; i++) {

        if (verbose) { bar.update(); }

        auto *toplesets = new index_t[n_vertices];
        auto *sorted_index = new index_t[n_vertices];
        vector<index_t> limits;
        vector<index_t> source = {static_cast<unsigned int>(i)};

        mesh->compute_toplesets(toplesets, sorted_index, limits, source);
        res[i] = fast_marching_single_vert(mesh, source);

        delete[] toplesets;
        delete[] sorted_index;

    }

    delete mesh;
    return res;
}



vector<double> fast_marching_single_vert(che *mesh, const vector<index_t> &source) {
    geodesics fm(mesh, source);
    vector<double> res;
    for (index_t v = 0; v < mesh->n_vertices(); v++) {
        res.push_back(fm[v]);
    }
    return res;
}

