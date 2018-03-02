#ifndef FFNN_WAVE_FUNCTION
#define FFNN_WAVE_FUNCTION

#include "WaveFunction.hpp"
#include "FeedForwardNeuralNetwork.hpp"

#include <stdexcept>



class FFNNWaveFunction: public WaveFunction{

private:
    FeedForwardNeuralNetwork * _ffnn;

    double _wf_value;
    double * _d1_logwf;
    double * _d2_logwf;
    double * _vd1_logwf;

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

    // --- computation of the derivatives
    void computeAllInternalValues(const double *in);

    // --- wf derivatives and other internal values
    double getD1LogWF(const int &id1){return _d1_logwf[id1];}
    double getD2LogWF(const int &id2){return _d2_logwf[id2];}
    double getVD1LogWF(const int &ivd1){return _vd1_logwf[ivd1];}
    double getWFValue(){return _wf_value;}

};


#endif
