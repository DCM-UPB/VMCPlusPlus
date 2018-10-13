#ifndef PAIR_SYMMETRIZER_WAVE_FUNCTION
#define PAIR_SYMMETRIZER_WAVE_FUNCTION


#include "SymmetrizerWaveFunction.hpp"


class PairSymmetrizerWaveFunction: public SymmetrizerWaveFunction {
    /*
      This class wraps around an arbitrary WaveFunction and applies by default the
      pair Symmetrizer operator or by optional second constructor argument instead
      the pair AntiSymmetrizer operator to the given wavefunction. These operators
      are like their general version, but use only permutations that arise from
      single transpositions of particle pairs. The result is a wavefunction that is
      only approximately symmetric or antisymmetric with respect to particle exchange.

      NOTE: Compared to the full Symmetrizer operation (with N! evaluations), this
      version uses only N*(N-1)/2+1 evaluations. However it is not necessarily fully
      symmetric/antisymmetric for permutations that involve more than one transposition.
    */
public:
    PairSymmetrizerWaveFunction(WaveFunction * wf, const bool flag_antisymmetric = false):
        SymmetrizerWaveFunction(wf, flag_antisymmetric){}

    ~PairSymmetrizerWaveFunction(){}

    // devirtualize inherited methods
    double getAcceptance(const double * protoold, const double * protonew) {
        return SymmetrizerWaveFunction::getAcceptance(protoold, protonew);
    }
    double computeWFValue(const double * protovalues) {
        return SymmetrizerWaveFunction::computeWFValue(protovalues);
    }

    // own implementations
    void samplingFunction(const double * in, double * out);
    void computeAllDerivatives(const double * x);
};

#endif
