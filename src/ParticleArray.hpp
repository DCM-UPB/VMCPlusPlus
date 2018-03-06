#ifndef PARTICLE_ARRAY
#define PARTICLE_ARRAY


class ParticleArray{

private:
    int _nspacedim;

public:
    ParticleArray(const int &nspacedim){
        _nspacedim = nspacedim;
    }

    double * getParticleArray(double * x, const int &i){return x+_nspacedim*i;}
    const double * getParticleArray(const double * x, const int &i){return x+_nspacedim*i;}

    void setParticleArray(double * x, const int &i, const double * newx){
        const int istart = _nspacedim*i;
        for (int k=0; k<_nspacedim; ++k){
            x[istart + k] = newx[k];
        }
    }

    void addArrayToParticleArray(double * x, const int &i, const double * toadd){
        const int istart = _nspacedim*i;
        for (int k=0; k<_nspacedim; ++k){
            x[istart + k] += toadd[k];
        }
    }

};


#endif
