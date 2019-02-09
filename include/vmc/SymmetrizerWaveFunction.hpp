#ifndef SYMMETRIZER_WAVE_FUNCTION
#define SYMMETRIZER_WAVE_FUNCTION


#include "vmc/WaveFunction.hpp"


class SymmetrizerWaveFunction: public WaveFunction {
    /*
      This class wraps around an arbitrary WaveFunction and applies by default the
      general Symmetrizer operator or by optional second constructor argument instead
      the general AntiSymmetrizer operator to the given wavefunction. The result is
      a wavefunction that is either symmetric or antisymmetric with respect to arbitrary
      particle exchange permutations.

      NOTE 1: These general Symmetrizer operators are not Slater determinants, so they don't
      rely on N single particle orbitals, but instead a single N-particle wavefunction.

      NOTE 2: Be aware that the time to evaluate the samplingFunction or derivatives of the
      SymmetrizerWaveFunction requires at least N! as much time as the corresponding evaluations
      of the original WaveFunction, where n is the number of particles. Therefore it cannot be used
      in practice for more than a handful of particles.
    */
protected:
    WaveFunction * _wf; // we wrap around an existing wavefunction
    const bool _flag_antisymmetric; // should we use the antisymmetrizer instead of symmetrizer?

    // internal helpers
    unsigned long _npart_factorial();
    void _swapPositions(double * x, const int &i, const int &j);
    void _swapIndices(int * ids, const int &i, const int &j);
    void _computeStandardDerivatives(const double * x, const double &normf);
    void _addSwapDerivatives(const double * x, const double &normf, const int * ids);

public:
    SymmetrizerWaveFunction(WaveFunction * wf, const bool flag_antisymmetric = false):
        WaveFunction(wf->getNSpaceDim(), wf->getNPart(), 1, wf->getNVP(), wf->hasVD1(), wf->hasD1VD1(), wf->hasD2VD1()), _wf(wf), _flag_antisymmetric(flag_antisymmetric) {}

    virtual ~SymmetrizerWaveFunction(){}

    void setVP(const double * vp);

    void getVP(double * vp);

    virtual void samplingFunction(const double * in, double * out);

    virtual double getAcceptance(const double * protoold, const double * protonew);

    virtual void computeAllDerivatives(const double * x);

    virtual double computeWFValue(const double * protovalues);

    bool isAntiSymmetric() {return _flag_antisymmetric;}
};


#endif
