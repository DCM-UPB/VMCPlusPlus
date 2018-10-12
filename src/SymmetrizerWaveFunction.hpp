#ifndef SYMMETRIZER_WAVE_FUNCTION
#define SYMMETRIZER_WAVE_FUNCTION


#include "WaveFunction.hpp"


class SymmetrizerWaveFunction: virtual public WaveFunction {
    /*
      This class wraps around an arbitrary WaveFunction and applies by default the
      general Symmetrizer operator or by optional second constructor argument instead
      the general AntiSymmetrizer operator to the given wavefunction. The result is
      a wavefunction that is either symmetric or antisymmetric with respect to particle
      exchange.

      NOTE 1: These general Symmetrizer operators are not Slater determinants, so they don't
      rely on single particle orbitals, but instead a single arbitrary n-particle wavefunction.

      NOTE 2: Be aware that the time to evaluate the samplingFunction or derivatives of the
      SymmetrizerWaveFunction requires at least n! as much time as the corresponding evaluations
      of the original WaveFunction, where n is the number of particles. Therefore it cannot be used
      in practice for more than a handful of particles.
    */
protected:
    WaveFunction * _wf; // we wrap around an existing wavefunction
    const bool _flag_antisymmetric; // should we use the antisymmetrizer instead of symmetrizer?

    // internal helpers
    size_t _npart_factorial();
    void _swapParticles(double * x, const int &i, const int &j);
    void _computeStandardDerivatives(const double * x, const double &normf);
    void _addSwapDerivatives(const double * x, const double &normf);

public:
    SymmetrizerWaveFunction(WaveFunction * wf, const bool flag_antisymmetric = false):
        WaveFunction(wf->getNSpaceDim(), wf->getNPart(), 1, wf->getNVP(), wf->hasVD1(), wf->hasD1VD1(), wf->hasD2VD1()), _wf(wf), _flag_antisymmetric(flag_antisymmetric) {}

    virtual ~SymmetrizerWaveFunction(){}

    void setVP(const double * vp);

    void getVP(double * vp);

    void samplingFunction(const double * in, double * out);

    double getAcceptance(const double * protoold, const double * protonew);

    void computeAllDerivatives(const double * x);

    double computeWFValue(const double * protovalues);

    bool isAntiSymmetric() {return _flag_antisymmetric;}
};


#endif
