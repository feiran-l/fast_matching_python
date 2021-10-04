#ifndef TEST_GEODESICS_PTP_H
#define TEST_GEODESICS_PTP_H

#include "geodesics.h"


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

vector<vector<double> > calc_geodesic_matrix(const string data_path, bool verbose=true);
vector<double> fast_marching_single_vert(che *mesh, const std::vector<index_t> &source);


#endif // TEST_GEODESICS_PTP_H

