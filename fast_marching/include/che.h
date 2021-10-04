#ifndef CHE_H
#define CHE_H

#include "include.h"
#include "vertex.h"


#define for_star(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->ot(prev(he))) != stop ? he : NIL)
#define for_border(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->evt(mesh->vt(next(he)))) != stop ? he : NIL)


typedef std::vector<index_t> link_t;        // link (vector of he)
index_t trig(const index_t &he);
index_t next(const index_t &he);
index_t prev(const index_t &he);


class che {
    public:
        static const size_t P = 3;
    protected:
        std::string filename_;
        size_t n_vertices_{};
        size_t n_faces_{};
        size_t n_half_edges_{};
        size_t n_edges_{};
        size_t n_borders_{};
        vertex *GT{};    ///< geometry table			: v		-> vertex
        index_t *VT{};    ///< vertex table (faces)	: he	-> v
        index_t *OT{};    ///< opposite table			: he	-> he
        index_t *EVT{};    ///< extra vertex table		: v		-> he
        index_t *ET{};    ///< edge table				: e		-> he
        index_t *EHT{};    ///< extra half edge table	: he	-> e
        index_t *BT{};    ///< boundary table			: b 	-> v
        bool manifold{};

public:
        che(const size_t &n_v = 0, const size_t &n_f = 0);
        che(const vertex *vertices, const index_t &n_v, const index_t *faces, const index_t &n_f);
        che(const che &mesh);
        virtual ~che();

        void link(link_t &l, const index_t &v);
        const index_t &vt(const index_t &he) const;
        const vertex &gt(const index_t &v) const;
        const index_t &ot(const index_t &he) const;
        const index_t &evt(const index_t &v) const;
        const size_t &n_vertices() const;
        const size_t &n_faces() const;
        const size_t &n_half_edges() const;
        void compute_toplesets(index_t *&rings, index_t *&sorted, std::vector<index_t> &limites,
                                const std::vector<index_t> &sources, const index_t &k = NIL);

    protected:
        void delete_me();
        void init(const vertex *vertices, const index_t &n_v, const index_t *faces, const index_t &n_f);
        void init(const std::string &file);
        void init(const size_t &n_v, const size_t &n_f);
        virtual void read_file(const std::string &file);

    private:
        void update_evt_ot_et();
        void update_eht();
        void update_bt();
    };




#endif // CHE_H

