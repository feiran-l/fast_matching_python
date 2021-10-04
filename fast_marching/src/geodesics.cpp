#include "geodesics.h"



geodesics::geodesics(che *mesh, const vector<index_t> &sources, double *const &e_dist, const bool &cluster, 
                        const size_t &n_iter, const double &radio) : n_vertices(mesh->n_vertices()) {
    assert(n_vertices > 0);
    free_dist = e_dist == nullptr;
    dist = free_dist ? new double[n_vertices] : e_dist;
    clusters = cluster ? new index_t[n_vertices] : nullptr;
    sorted_index = new index_t[n_vertices];
    n_sorted = 0;

    memset(sorted_index, -1, n_vertices * sizeof(index_t));
    for (index_t v = 0; v < n_vertices; v++)
        dist[v] = INFINITY;

    assert(!sources.empty());
    run_fast_marching(mesh, sources, n_iter, radio);
}


geodesics::~geodesics() {
    if (free_dist) delete[] dist;
    delete[] sorted_index;
    delete[] clusters;
}


const double &geodesics::operator[](const index_t &i) const {
    assert(i < n_vertices);
    return dist[i];
}


const index_t &geodesics::operator()(const index_t &i) const {
    assert(i < n_vertices);
    return sorted_index[i];
}


const index_t &geodesics::farthest() const {
    assert(n_sorted != 0);
    return sorted_index[n_sorted - 1];
}




void geodesics::run_fast_marching(che *mesh, const vector<index_t> &sources, const size_t &n_iter, const double &radio) {
    index_t BLACK = 0, GREEN = 1, RED = 2;
    auto *color = new index_t[n_vertices];

    for (index_t v = 0; v < n_vertices; v++)
        color[v] = GREEN;

    size_t green_count = n_iter ? n_iter : n_vertices;
    priority_queue<pair<double, size_t>, vector<pair<double, size_t> >, greater<pair<double, size_t> > > Q;

    double dv, dp;
    index_t dir; // dir propagation
    vertex vx;
    size_t black_i, v;
    index_t c = 0;
    n_sorted = 0;

    for (index_t s: sources) {
        dist[s] = 0;
        if (clusters) clusters[s] = ++c;
        color[s] = RED;
        Q.push(make_pair(dist[s], s));
    }

    while (green_count-- && !Q.empty()) {
        while (!Q.empty() && color[Q.top().second] == BLACK)
            Q.pop();

        if (Q.empty()) break;

        black_i = Q.top().second;
        color[black_i] = BLACK;
        Q.pop();

        if (dist[black_i] > radio) break;

        sorted_index[n_sorted++] = black_i;

        link_t black_link;
        mesh->link(black_link, black_i);
        for (const index_t &he: black_link) {
            v = mesh->vt(he);

            if (color[v] == GREEN)
                color[v] = RED;

            if (color[v] == RED) {
                dv = dist[v];
                for_star(v_he, mesh, v) {
                    dp = update(dir, mesh, v_he, vx);
                    if (dp < dv) {
                        dv = dp;

                        if (clusters)
                            clusters[v] = dist[mesh->vt(prev(v_he))] < dist[mesh->vt(next(he))] ? clusters[mesh->vt(
                                    prev(he))] : clusters[mesh->vt(next(he))];
                    }
                }

                if (dv < dist[v])
                    Q.push(make_pair(dist[v] = dv, v));
            }
        }
    }

    delete[] color;
}





double geodesics::update(index_t &d, che *mesh, const index_t &he, vertex &vx) {
    d = NIL;

    Matrix<double, 3, 2> X;

    index_t x[3];

    x[0] = mesh->vt(next(he));
    x[1] = mesh->vt(prev(he));
    x[2] = mesh->vt(he);                //update x[2]

    vx = mesh->gt(x[2]);

    vertex v[2];
    v[0] = mesh->gt(x[0]) - vx;
    v[1] = mesh->gt(x[1]) - vx;

    X(0, 0) = v[0][0]; X(1, 0) = v[0][1]; X(2, 0) = v[0][2];
    X(0, 1) = v[1][0]; X(1, 1) = v[1][1]; X(2, 1) = v[1][2];

    return planar_update(d, X, x, vx);
}




double geodesics::planar_update(index_t &d, Matrix<double, 3, 2> &X, index_t *x, vertex &vx) {

    Vector2d ones {1, 1};

    Matrix2d Q;
    LLT<Matrix2d> llt_Q(X.transpose() * X); // compute the Cholesky decomposition of A
    if(llt_Q.info() == Eigen::NumericalIssue) 
        return INFINITY;
    else
        Q = (X.transpose() * X).inverse();


    Vector2d t {dist[x[0]], dist[x[1]]};

    double p;
    double delta = ones.transpose() * Q * t;

    double scalar1 = ones.transpose() * Q * ones;
    double scalar2 = t.transpose() * Q * t - 1;
    double dis = delta * delta - scalar1 * scalar2 ;
    // distance_t dis = delta * delta - (ones.transpose() * Q * ones) * (t.transpose() * Q * t - 1);

    if (dis >= 0) {
        p = delta + sqrt(dis);
        p /= ones.transpose() * Q * ones;
    } else p = INFINITY;

    MatrixXd n = X * Q * (t - p * ones);
    MatrixXd cond = Q * X.transpose() * n;
    Vector3d v(3);

    if (t(0) == INFINITY || t(1) == INFINITY || dis < 0 || (cond(0) >= 0 || cond(1) >= 0)) {
        double dp[2];
        dp[0] = dist[x[0]] + X.col(0).norm();
        dp[1] = dist[x[1]] + X.col(1).norm();
        d = dp[1] < dp[0];
        v = X.col(d);
        p = dp[d];
    } 
    else {
        Matrix<double, 3, 2> A;
        A.col(0) = -n;
        A.col(1) = X.col(1) - X.col(0);
        Vector3d b = -X.col(0);
        MatrixXd l = A.colPivHouseholderQr().solve(b);
        v = l(1) * A.col(1) + X.col(0);
    }

    vx += vertex(v(0), v(1), v(2));

    return p;
}


