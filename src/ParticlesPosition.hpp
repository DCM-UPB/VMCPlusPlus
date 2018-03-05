#ifndef PARTICLES_POSITION
#define PARTICLES_POSITION


class ParticlesPosition{

private:
    int _nspacedim;

public:
    ParticlePosition(const int &nspacedim){
        _nspacedim = nspacedim;
    }

    double * getParticlePosition(const double *x, const int &i){return x+_nspacedim*i;}

    void setParticlePosition(double * x, const int &i, const double * newx){
        const int istart = _nsapcedim*i;
        for (int k=0; k<_nsapcedim; ++k){
            x[istart + k] = newx[k];
        }
    }

    void addArrayToParticlePosition(double * x, const int &i, const double * toadd){
        const int istart = _nsapcedim*i;
        for (int k=0; k<_nsapcedim; ++k){
            x[istart + k] += newx[k];
        }
    }

};


#endif
