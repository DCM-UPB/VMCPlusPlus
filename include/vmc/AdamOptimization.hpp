#ifndef VMC_ADAMOPTIMIZATION_HPP
#define VMC_ADAMOPTIMIZATION_HPP

#include "nfm/Adam.hpp"
#include "vmc/EnergyGradientTargetFunction.hpp"
#include "vmc/StochasticReconfigurationTargetFunction.hpp"
#include "vmc/WFOptimization.hpp"

#include "mci/MCIntegrator.hpp"

namespace vmc
{

class AdamOptimization: public WFOptimization
{
private:
    const int _Nmc;
    const bool _useSR;
    const bool _useGradientError;
    const bool _useAveraging;
    const int _max_n_const_values;
    const double _lambda;
    const double _alpha;
    const double _beta1, _beta2;
    const double _epsilon;

public:
    AdamOptimization(WaveFunction * wf, Hamiltonian * H, mci::MCI * mci, int Nmc, bool useSR = false, bool useGradientError = false, int max_n_const_values = 20,
                     bool useAveraging = false, double lambda = 0., double alpha = 0.001, double beta1 = 0.9, double beta2 = 0.999, double epsilon = 10e-8):
            WFOptimization(wf, H, mci),
            _Nmc(Nmc), _useSR(useSR), _useGradientError(useGradientError), _useAveraging(useAveraging), _max_n_const_values(max_n_const_values),
            _lambda(lambda), _alpha(alpha), _beta1(beta1), _beta2(beta2), _epsilon(epsilon) {}

    ~AdamOptimization() override = default;

    // optimization
    void optimizeWF() override
    {
        // create targetfunction
        std::unique_ptr<nfm::NoisyFunctionWithGradient> targetf;
        if (_useSR) {
            targetf = std::make_unique<StochasticReconfigurationTargetFunction>(_wf, _H, getMCI(), _Nmc, _lambda, _useGradientError);
        }
        else { targetf = std::make_unique<EnergyGradientTargetFunction>(_wf, _H, _Nmc, _Nmc, getMCI(), _useGradientError, _lambda); }

        // declare the Adam object
        nfm::Adam adam(targetf->getNDim(), _useAveraging, _alpha);
        adam.setMaxNConstValues(_max_n_const_values);
        adam.setBeta1(_beta1);
        adam.setBeta2(_beta2);
        adam.setEpsilon(_epsilon);

        // allocate an array that will contain the wave function variational parameters
        double wfpar[_wf->getNVP()];
        // get the variational parameters
        _wf->getVP(wfpar);
        // set the actual variational parameters as starting point for optimization
        adam.setX(wfpar);
        // find the optimal parameters by minimizing the energy with the Adam algorithm
        adam.findMin(*targetf);
        // set the found parameters in the wave function
        _wf->setVP(adam.getX().data());
    }
};
} // namespace vmc


#endif
