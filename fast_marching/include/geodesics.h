#ifndef GEODESICS_H
#define GEODESICS_H

#include "che.h"
#include "include.h"



/*!
	Compute the geodesics distances on a mesh from a source or multi-source. This class implements
	the Fast Marching algorithm without deal with obtuse triangles. Also, if the options PTP_CPU or
	PTP_GPU are enables, compute the geodesics distances executing the Parallel Toplesets Propagation
	algorithm.
*/
class geodesics {

    public:
        index_t *clusters;            ///< Clustering vertices to closest source.

    private:
        double *dist;            ///< Results of computation geodesic distances.
        index_t *sorted_index;        ///< Sort vertices by topological level or geodesic distance.
        const size_t &n_vertices;    ///< Number of vertices.
        size_t n_sorted;            ///< Number of vertices sorted by their geodesics distance.
        bool free_dist;

    public:
        geodesics(che *mesh,                            ///< input mesh must be a triangular mesh.
                  const std::vector<index_t> &sources,    ///< source vertices.
                  double *const &e_dist = nullptr,    ///< external dist allocation
                  const bool &cluster = false,            ///< to cluster vertices to closest source.
                  const size_t &n_iter = 0,                ///< maximum number of iterations.
                  const double &radio = INFINITY        ///< execute until the specific radio.
        );

        virtual ~geodesics();

        const double &operator[](const index_t &i) const;

        const index_t &operator()(const index_t &i) const;


        const index_t &farthest() const;

    private:
        void run_fast_marching(che *mesh, const std::vector<index_t> &sources, const size_t &n_iter, const double &radio);
        double update(index_t &d, che *mesh, const index_t &he, vertex &vx);
        double planar_update(index_t &d, Matrix<double, 3, 2> &X, index_t *x, vertex &vx);
    };



#endif //GEODESICS_H

