#include "EuclidianParticleDistance.hpp"

#include <math.h>


double EuclidianParticleDistance::dist(const double * r1, const double * r2){
    double dist = pow(r1[0] - r2[0], 2);
    for (int i=1; i<getNSpaceDim(); ++i){
        dist += pow(r1[i] - r2[i], 2);
    }
    return sqrt(dist);
}

void EuclidianParticleDistance::distD1(const double * r1, const double * r2, double * out){
    const double invr = 1./dist(r1, r2);
    for (int i=0; i<getNSpaceDim(); ++i){
        out[i] = (r1[i] - r2[i]) * invr;
    }
}

void EuclidianParticleDistance::distD2(const double * r1, const double * r2, double * out){
    const double r = dist(r1, r2);
    for (int i=0; i<getNSpaceDim(); ++i){
        out[i] = ( pow(r,2) - pow(r1[i]-r2[i], 2) ) / pow(r, 3);
    }
}
