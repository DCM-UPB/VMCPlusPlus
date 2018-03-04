#ifndef HAMILTONIAN
#define HAMILTONIAN

#include "WaveFunction.hpp"
#include "MCIObservableFunctionInterface.hpp"


// Hamiltonian: obs[0]=Totale Energy, obs[1]=Potential Energy, obs[2]=Kinetic Energy (PB), obs[3]=Kinetic Energy (JF)
class Hamiltonian: public MCIObservableFunctionInterface
{
protected:
    int _nspacedim;
    int _npart;
    WaveFunction * _wf;

public:
    Hamiltonian(const int &nspacedim, const int &npart, WaveFunction * wf): MCIObservableFunctionInterface(nspacedim*npart, 4)
    {
        _nspacedim=nspacedim; _npart=npart; _wf = wf;
    }
    virtual ~Hamiltonian(){}

    int getNSpaceDim(){return _nspacedim;}
    int getTotalNDim(){return getNDim();}
    int getNPart(){return _npart;}

    // Potential energy --- MUST BE IMPLEMENTED
    virtual double localPotentialEnergy(const double *r) = 0;


    double localPBKineticEnergy(const double *r)
    {
        double ekin=0.;
        for (int i=0; i<_ndim; ++i)
            {
                ekin += _wf->getD2LogWF(i);
            }
        return (-0.5*ekin);
    }

    double localJFKineticEnergy(const double *r)
    {
        double ekin=0.;
        double foo;
        for (int i=0; i<_ndim; ++i)
            {
                foo = _wf->getD1LogWF(i);
                ekin += foo*foo;
            }
        return (0.5*ekin);
    }

    void observableFunction(const double * in, double *out)
    {
        out[3]=localJFKineticEnergy(in);
        out[2]=localPBKineticEnergy(in);
        out[1]=localPotentialEnergy(in);
        out[0]=out[1]+out[2];
    }
};


#endif
