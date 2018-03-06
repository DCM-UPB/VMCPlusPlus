#ifndef EUCLIDIAN_PARTICLE_DISTANCE
#define EUCLIDIAN_PARTICLE_DISTANCE

#include "ParticlesDistance.hpp"


class EuclidianParticlesDistance: public ParticlesDistance{

public:
    EuclidianParticlesDistance(const int &nspacedim): ParticlesDistance(nspacedim){}


    double dist(const double * r1, const double * r2);

    void distD1(const double * r1, const double * r2, double * out);

    void distD2(const double * r1, const double * r2, double * out);

};


#endif
