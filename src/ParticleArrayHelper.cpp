#include "vmc/ParticleArrayHelper.hpp"

#include <stdexcept>



double * ParticleArrayHelper::getParticleArray(double * x, const int &i){
    return x+_nspacedim*i;
}


const double * ParticleArrayHelper::getParticleArray(const double * x, const int &i){
    return x+_nspacedim*i;
}


void ParticleArrayHelper::setParticleArray(double * x, const int &i, const double * newx){
    const int istart = _nspacedim*i;
    for (int k=0; k<_nspacedim; ++k){
        x[istart + k] = newx[k];
    }
}


void ParticleArrayHelper::addArrayToParticleArray(double * x, const int &i, const double * toadd){
    const int istart = _nspacedim*i;
    for (int k=0; k<_nspacedim; ++k){
        x[istart + k] += toadd[k];
    }
}
