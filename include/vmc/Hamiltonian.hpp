#ifndef VMC_HAMILTONIAN_HPP
#define VMC_HAMILTONIAN_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "vmc/WaveFunction.hpp"

namespace vmc
{

// Hamiltonian: obs[0]=Totale Energy, obs[1]=Potential Energy, obs[2]=Kinetic Energy (PB), obs[3]=Kinetic Energy (JF)
class Hamiltonian: public mci::ObservableFunctionInterface
{
protected:
    const int _nspacedim;
    const int _npart;
    WaveFunction * const _wf;

    const bool _flag_PBKE;

public:
    Hamiltonian(const int &nspacedim, const int &npart, WaveFunction * wf, const bool usePBKE = true /* use JF+PB KE or only JF */)
        : mci::ObservableFunctionInterface(nspacedim*npart, 4), _nspacedim(nspacedim), _npart(npart), _wf(wf), _flag_PBKE(usePBKE) {}
    ~Hamiltonian() override= default;

    int getNSpaceDim(){return _nspacedim;}
    int getTotalNDim(){return getNDim();}
    int getNPart(){return _npart;}

    bool hasPBKE(){return _flag_PBKE;}

    // Potential energy --- MUST BE IMPLEMENTED
    virtual double localPotentialEnergy(const double *r) = 0;


    double localPBKineticEnergy(const double * /*r*/)
    {
        double ekin=0.;
        for (int i=0; i<_ndim; ++i)
            {
                ekin += _wf->getD2DivByWF(i);
            }
        return (-0.5*ekin);
    }

    double localJFKineticEnergy(const double * /*r*/)
    {
        double ekin=0.;
        for (int i=0; i<_ndim; ++i)
            {
                const double foo = _wf->getD1DivByWF(i);
                ekin += foo*foo;
            }
        return (0.5*ekin);
    }

    void observableFunction(const double * in, double *out) override
    {
        out[3]=localJFKineticEnergy(in);
        out[2]= _flag_PBKE ? localPBKineticEnergy(in) : 0.;
        out[1]=localPotentialEnergy(in);
        out[0]= _flag_PBKE ? out[1]+out[2] : out[1]+out[3];
    }
};
} // namespace vmc

#endif
