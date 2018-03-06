#include "ParticleArrayManager.hpp"

#include <stdexcept>



double * ParticleArrayManager::getParticleArray(double * x, const int &i){
    return x+_nspacedim*i;
}


const double * ParticleArrayManager::getParticleArray(const double * x, const int &i){
    return x+_nspacedim*i;
}


void ParticleArrayManager::setParticleArray(double * x, const int &i, const double * newx){
    const int istart = _nspacedim*i;
    for (int k=0; k<_nspacedim; ++k){
        x[istart + k] = newx[k];
    }
}


void ParticleArrayManager::addArrayToParticleArray(double * x, const int &i, const double * toadd){
    const int istart = _nspacedim*i;
    for (int k=0; k<_nspacedim; ++k){
        x[istart + k] += toadd[k];
    }
}
