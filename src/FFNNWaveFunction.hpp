#ifndef FFNN_WAVE_FUNCTION
#define FFNN_WAVE_FUNCTION

#include "WaveFunction.hpp"
#include "FeedForwardNeuralNetwork.hpp"

#include <stdexcept>



class FFNNWaveFunction: public WaveFunction{

private:
    FeedForwardNeuralNetwork * _ffnn;

public:
    // --- Constructor and destructor
    // IMPORTANT: The provided ffnn should be ready to use (connected) and have the first, second and variational derivatives substrates
    FFNNWaveFunction(const int &nspacedim, const int &npart, FeedForwardNeuralNetwork * ffnn);
    ~FFNNWaveFunction();



    // --- Getters
    FeedForwardNeuralNetwork * getFFNN(){return _ffnn;}


    // --- interface for manipulating the variational parameters
    void setVP(const double *vp);
    void getVP(double *vp);

    // --- methods herited from MCISamplingFunctionInterface
    // wave function values that will be used to compute the acceptance
    void samplingFunction(const double * in, double * out);
    // MCI acceptance starting from the new and old sampling functions
    double getAcceptance();

    // --- wf derivatives
    // first derivative divided by the wf
    double d1(const int &, const double * ); // 0 <= i <= _ndim
    // second derivative divided by the wf
    double d2(const int &, const double *);  // 0 <= i <= _ndim  ,  0 <= j <= _ndim
    // variational derivative divided by the wf
    double vd1(const int &, const double *);  // 0 <= i <= _npv
};


#endif
