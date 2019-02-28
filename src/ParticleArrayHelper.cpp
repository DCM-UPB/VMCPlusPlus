#include "vmc/ParticleArrayHelper.hpp"

#include <algorithm>



double * ParticleArrayHelper::getParticleArray(double * x, const int &i){
    return x+i*_nspacedim;
}


const double * ParticleArrayHelper::getParticleArray(const double * x, const int &i){
    return x+i*_nspacedim;
}


void ParticleArrayHelper::setParticleArray(double * x, const int &i, const double * newx){
    std::copy(newx, newx+_nspacedim, x+i*_nspacedim);
}


void ParticleArrayHelper::addArrayToParticleArray(double * x, const int &i, const double * toadd){
    const int istart = i*_nspacedim;
    for (int k=0; k<_nspacedim; ++k){
        x[istart + k] += toadd[k];
    }
}
