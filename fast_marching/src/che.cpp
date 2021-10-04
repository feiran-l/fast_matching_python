#include "che.h"



index_t trig(const index_t &he) {
    if (he == NIL) return NIL;
    return he / che::P;
}

index_t next(const index_t &he) {
    if (he == NIL) return NIL;
    return che::P * trig(he) + (he + 1) % che::P;
}

index_t prev(const index_t &he) {
    if (he == NIL) return NIL;
    return che::P * trig(he) + (he + che::P - 1) % che::P;
}


/*-----------------------------------------------------------------------------------------------------------*/


che::che(const che &mesh) {
    filename_ = mesh.filename_;
    n_vertices_ = mesh.n_vertices_;
    n_faces_ = mesh.n_faces_;
    n_half_edges_ = mesh.n_half_edges_;
    n_edges_ = mesh.n_edges_;
    n_borders_ = mesh.n_borders_;

    GT = new vertex[n_vertices_];
    memcpy(GT, mesh.GT, n_vertices_ * sizeof(vertex));

    VT = new index_t[n_half_edges_];
    memcpy(VT, mesh.VT, n_half_edges_ * sizeof(index_t));

    OT = new index_t[n_half_edges_];
    memcpy(OT, mesh.OT, n_half_edges_ * sizeof(index_t));

    EVT = new index_t[n_vertices_];
    memcpy(EVT, mesh.EVT, n_vertices_ * sizeof(index_t));

    ET = new index_t[n_edges_];
    memcpy(ET, mesh.ET, n_edges_ * sizeof(index_t));

    EHT = new index_t[n_half_edges_];
    memcpy(EHT, mesh.EHT, n_half_edges_ * sizeof(index_t));

    BT = new index_t[n_borders_];
    memcpy(BT, mesh.BT, n_borders_ * sizeof(index_t));
}

che::che(const size_t &n_v, const size_t &n_f) {
    init(n_v, n_f);
}

che::che(const vertex *vertices, const index_t &n_v, const index_t *faces, const index_t &n_f) {
    init(vertices, n_v, faces, n_f);
}

che::~che() {
    delete_me();
}


void che::link(link_t &l, const index_t &v) {
    if (v >= n_vertices_) return;
    for_star(he, this, v) {
        l.push_back(next(he));
        if (OT[prev(he)] == NIL)
            l.push_back(prev(he));
    }
}


const index_t &che::vt(const index_t &he) const {
    assert(he < n_half_edges_);
    return VT[he];
}


const vertex &che::gt(const index_t &v) const {
    assert(v < n_vertices_);
    return GT[v];
}


const index_t &che::ot(const index_t &he) const {
    assert(he < n_half_edges_);
    return OT[he];
}

const index_t &che::evt(const index_t &v) const {
    assert(v < n_vertices_);
    return EVT[v];
}

const size_t &che::n_vertices() const {
    return n_vertices_;
}

const size_t &che::n_faces() const {
    return n_faces_;
}

const size_t &che::n_half_edges() const {
    return n_half_edges_;
}


void che::compute_toplesets(index_t *&toplesets, index_t *&sorted, vector<index_t> &limits,
                            const vector<index_t> &sources, const index_t &k) {
    if (sources.empty()) return;

    memset(toplesets, -1, sizeof(index_t) * n_vertices_);

    index_t level = 0;

    index_t p = 0;
    for (const index_t &s: sources) {
        sorted[p++] = s;
        if (toplesets[s] == NIL)
            toplesets[s] = level;
    }

    limits.push_back(0);
    for (index_t i = 0; i < p; i++) {
        const index_t &v = sorted[i];
        if (toplesets[v] > level) {
            level++;
            if (level > k) break;
            limits.push_back(i);
        }

        link_t v_link;
        link(v_link, v);
        for (const index_t &he: v_link) {
            const index_t &u = VT[he];
            if (toplesets[u] == NIL) {
                toplesets[u] = toplesets[v] + 1;
                sorted[p++] = u;
            }
        }
    }

    assert(p <= n_vertices_);
    limits.push_back(p);
}


void che::delete_me() {
    delete[] GT;
    delete[] VT;
    delete[] OT;
    delete[] EVT;
    delete[] ET;
    delete[] EHT;
    delete[] BT;
}


void che::init(const vertex *vertices, const index_t &n_v, const index_t *faces, const index_t &n_f) {
    init(n_v, n_f);
    memcpy(GT, vertices, n_vertices_ * sizeof(vertex));
    memcpy(VT, faces, n_half_edges_ * sizeof(index_t));
    update_evt_ot_et();
    update_eht();
    update_bt();
}


void che::init(const string &file) {
    filename_ = file;
    read_file(filename_);
    update_evt_ot_et();
    update_eht();
    update_bt();
}


void che::init(const size_t &n_v, const size_t &n_f) {
    n_vertices_ = n_v;
    n_faces_ = n_f;
    n_half_edges_ = n_edges_ = n_borders_ = 0;
    GT = nullptr;
    VT = OT = EVT = ET = BT = nullptr;
    manifold = true;
    n_half_edges_ = che::P * n_faces_;
    n_edges_ = 0; //n_half_edges_ / 2;	/**/
    if (n_vertices_) GT = new vertex[n_vertices_];
    if (n_half_edges_) VT = new index_t[n_half_edges_];
    if (n_half_edges_) OT = new index_t[n_half_edges_];
    if (n_vertices_) EVT = new index_t[n_vertices_];
    if (n_vertices_) EHT = new index_t[n_half_edges_];
}


void che::update_evt_ot_et() {
    memset(EVT, -1, sizeof(index_t) * n_vertices_);
    if (!n_faces_) return;

    auto *he_p_vertex = new vector<index_t>[n_vertices_];

    //vertex table
    for (index_t he = 0; he < n_half_edges_; he++) {
        EVT[VT[he]] = he;
        he_p_vertex[VT[he]].push_back(he);
    }

    //opposite table - edge table
    memset(OT, -1, sizeof(index_t) * n_half_edges_);

    vector<index_t> et;
    for (index_t he = 0; he < n_half_edges_; he++) {
        if (OT[he] == NIL) {
            et.push_back(he);
            for (index_t h: he_p_vertex[VT[he]]) {
                if (VT[prev(h)] == VT[next(he)]) {
                    if (OT[he] == NIL && OT[prev(h)] == NIL) {
                        OT[he] = prev(h);
                        OT[prev(h)] = he;
                    }
                };
            };
        }
    }

    //edge table
    n_edges_ = et.size();
    ET = new index_t[n_edges_];
    memcpy(ET, et.data(), sizeof(index_t) * n_edges_);

    for (index_t he = 0; he < n_half_edges_; he++)
        if (OT[he] == NIL && EVT[VT[he]] != NIL) {
            if (OT[EVT[VT[he]]] == NIL && EVT[VT[he]] != he) {
                manifold = false;
                EVT[VT[he]] = NIL;
            } else EVT[VT[he]] = he;
        }

    delete[] he_p_vertex;
}

void che::update_eht() {
    for (index_t e = 0; e < n_edges_; e++) {
        EHT[ET[e]] = e;
        if (OT[ET[e]] != NIL)
            EHT[OT[ET[e]]] = e;
    }
}


void che::update_bt() {
    if (!n_faces_) return;
    if (!manifold) return;

    bool *border = new bool[n_vertices_];
    memset(border, 0, sizeof(bool) * n_vertices_);

    vector<index_t> borders;

    for (index_t v = 0; v < n_vertices_; v++)
        if (!border[v] && EVT[v] != NIL && OT[EVT[v]] == NIL) {
            borders.push_back(v);
            for_border(he, this, v) border[VT[he]] = true;
        }

    n_borders_ = borders.size();
    if (n_borders_) {
        BT = new index_t[n_borders_];
        memcpy(BT, borders.data(), sizeof(index_t) * n_borders_);
    } else BT = nullptr;

    delete[] border;
}


void che::read_file(const string &) { }



