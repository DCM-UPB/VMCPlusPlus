#ifndef EUCLIDIAN_PARTICLE_DISTANCE
#define EUCLIDIAN_PARTICLE_DISTANCE

#include "ParticleDistance.hpp"


class EuclidianParticleDistance: public ParticleDistance{

public:
    EuclidianParticleDistance(const int &nspacedim): ParticleDistance(nspacedim){}
    ~EuclidianParticleDistance(){}

    double dist(const double * r1, const double * r2);

    void distD1(const double * r1, const double * r2, double * out);

    void distD2(const double * r1, const double * r2, double * out);

};


#endif
