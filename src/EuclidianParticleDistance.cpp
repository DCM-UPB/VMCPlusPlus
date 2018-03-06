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
        out[getNSpaceDim()+i] = (r2[i] - r1[i]) * invr;
    }
}

void EuclidianParticleDistance::distD2(const double * r1, const double * r2, double * out){
    const double r = dist(r1, r2);
    const double rpow2 = pow(r, 2);
    const double rpow3 = r*rpow2;
    for (int i=0; i<getNSpaceDim(); ++i){
        const double ripow2 = pow(r1[i]-r2[i], 2);
        out[i] = ( rpow2 - ripow2 ) / rpow3;
        out[getNSpaceDim()+i] = ( rpow2 - ripow2 ) / rpow3;
    }
}
