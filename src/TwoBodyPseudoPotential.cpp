#include "TwoBodyPseudoPotential.hpp"



double TwoBodyPseudoPotential::u(const double * r1, const double * r2){
    return ur(_dist->dist(r1, r2));
}


void TwoBodyPseudoPotential::d1(const double * r1, const double * r2, double * deriv1){
    const double ud1 = urD1(_dist->dist(r1, r2));
    _dist->distD1(r1, r2, deriv1);
    for (int i=0; i<2*_dist->getNSpaceDim(); ++i){
        deriv1[i] *= ud1;
    }
}


void TwoBodyPseudoPotential::d2(const double * r1, const double * r2, double * deriv2){
    const double ud1 = urD1(_dist->dist(r1, r2));
    const double ud2 = urD2(_dist->dist(r1, r2));

    _dist->distD2(r1, r2, deriv2);
    _dist->distD1(r1, r2, _foo);

    for (int i=0; i<2*_dist->getNSpaceDim(); ++i){
        deriv2[i] *= deriv2[i] * ud1 + _foo[i] * _foo[i] * ud2;
    }
}
