#include "vmc/EuclideanMetric.hpp"

#include <cmath>


double EuclideanMetric::dist(const double * r1, const double * r2){
    double dist = 0.;
    for (int i=0; i<getNSpaceDim(); ++i){
        dist += (r1[i] - r2[i])*(r1[i] - r2[i]);
    }
    return sqrt(dist);
}

void EuclideanMetric::distD1(const double * r1, const double * r2, double * out){
    const double invr = 1./dist(r1, r2);
    for (int i=0; i<getNSpaceDim(); ++i){
        out[i] = (r1[i] - r2[i]) * invr;
        out[getNSpaceDim()+i] = (r2[i] - r1[i]) * invr;
    }
}

void EuclideanMetric::distD2(const double * r1, const double * r2, double * out){
    const double r = dist(r1, r2);
    const double rpow2 = r*r;
    const double rpow3 = r*rpow2;
    for (int i=0; i<getNSpaceDim(); ++i){
        const double ripow2 = (r1[i]-r2[i])*(r1[i]-r2[i]);
        out[i] = ( rpow2 - ripow2 ) / rpow3;
        out[getNSpaceDim()+i] = ( rpow2 - ripow2 ) / rpow3;
    }
}
