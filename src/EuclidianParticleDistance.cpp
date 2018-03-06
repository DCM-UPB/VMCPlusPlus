#include "EuclidianParticleDistance.hpp"

#include <math.h>


double EuclidianParticleDistance::dist(const double * r1, const double * r2){
    double dist = pow(r2[0] - r1[0], 2);
    for (int i=1; i<getNSpaceDim(); ++i){
        dist += pow(r2[i] - r1[i], 2);
    }
    return sqrt(dist);
}

void EuclidianParticleDistance::distD1(const double * r1, const double * r2, double * out){
    const double invr = 1./dist(r1, r2);
    for (int i=0; i<getNSpaceDim(); ++i){
        out[i] = (r2[i] - r1[i]) * invr;
    }
}

void EuclidianParticleDistance::distD2(const double * r1, const double * r2, double * out){
    const double r = dist(r1, r2);
    const double invr5 = 1./pow(r, 5);
    const double r2times3 = 3.*pow(r, 2);
    for (int i=0; i<getNSpaceDim(); ++i){
        out[i] = (r2[i] - r1[i]) * invr5 * ( pow(r2[i] - r1[i], 2) - r2times3 );
    }
}
