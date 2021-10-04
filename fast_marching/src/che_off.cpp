#include "che_off.h"




che_off::che_off(const string &file) {
    init(file);
}


che_off::~che_off() {}



void che_off::read_file(const string &file) {
    string soff;
    size_t n_v, n_f, v;

    ifstream is(file);

    assert(is.good());

    is >> soff;
    is >> n_v >> n_f >> v;
    init(n_v, n_f);

    int r, g, b, a;
    for (index_t i = 0; i < n_vertices_; i++) {
        is >> GT[i];
        if (soff[0] == 'C') // COFF file, ignore RGBA
            is >> r >> g >> b >> a;
        if (soff[0] == 'N') // NOFF file, ignore normals
            is >> r >> g >> b;
    }

    index_t he = 0;
    for (index_t i = 0; i < n_faces_; i++) {
        is >> v;
        if (!i && v > che::P) {
            vertex *tGT = GT;
            GT = nullptr;
            delete_me();
            init(n_v, n_f * (v - che::P + 1));
            GT = tGT;
        }

        for (index_t j = 0; j < v; j++)
            is >> VT[he++];

        // divide face
        if (v == 4) {
            VT[he] = VT[he - v];
            he++;
            VT[he] = VT[he - che::P];
            he++;
            i++;
        }
    }

    is.close();
}




