#include <iostream>
#include <cmath>
#include <stdexcept>

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"
#include "ConjGrad.hpp"

#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"



/*
Hamiltonian describing a 1-particle harmonic oscillator:
   H  =  p^2 / 2m  +  1/2 * w^2 * x^2 
*/
class HarmonicOscillator1D1P: public Hamiltonian{

protected:
   double _w;

public:
   HarmonicOscillator1D1P(const double w, WaveFunction * wf): 
      Hamiltonian(1 /*num space dimensions*/, 1 /*num particles*/, wf) {_w=w;}
   
   // potential energy
   double localPotentialEnergy(const double *r)
   {
      return (0.5*_w*_w*(*r)*(*r));
   }
};



/*
Trial Wave Function made from a Feed Forward Neural Network
*/
class NeuralWaveFunction: public WaveFunction{
   private:
      static FeedForwardNeuralNetwork * buildNeuralNetwork(const int nhiddenlayers, const int hiddenlayersize){
         // check input parameters
         if (nhiddenlayers < 1){
            throw std::invalid_argument( "nhiddenlayers must be >= 1" );
         }
         if (hiddenlayersize < 2){
            throw std::invalid_argument( "hiddenlayersize must be >= 2" );
         }
         // generate FFNN
         FeedForwardNeuralNetwork * ffnn = new FeedForwardNeuralNetwork(2, hiddenlayersize, 2);
         for (int i=0; i<nhiddenlayers-1; ++i){
            ffnn->pushHiddenLayer(hiddenlayersize);
         }
         ffnn->connectFFNN();
         return ffnn;
      }
      
      static int countNeuralNetworkBeta(const int nhiddenlayers, const int hiddenlayersize){
         FeedForwardNeuralNetwork * ffnn = buildNeuralNetwork(nhiddenlayers, hiddenlayersize);
         const int result = ffnn->getNBeta();
         ffnn->disconnectFFNN();
         delete ffnn;
         return result;
      }
   
   protected:
      FeedForwardNeuralNetwork * _ffnn;

   public:
      NeuralWaveFunction(const int nhiddenlayers, const int hiddenlayersize): 
         WaveFunction(1 /*num space dimensions*/, 
            1 /*num particles*/, 
            1 /*num wf components*/, 
            countNeuralNetworkBeta(nhiddenlayers, hiddenlayersize) /*num variational parameters*/) {
               _ffnn = buildNeuralNetwork(nhiddenlayers, hiddenlayersize);
               _ffnn->addFirstDerivativeSubstrate();
               _ffnn->addSecondDerivativeSubstrate();
               _ffnn->addVariationalFirstDerivativeSubstrate();
            }
            
      ~NeuralWaveFunction(){
         _ffnn->disconnectFFNN();
         delete _ffnn;
      }
      
      
      FeedForwardNeuralNetwork * getFFNN(){return _ffnn;}
      
      
      void setVP(const double *in){
         for (int i=0; i<_ffnn->getNBeta(); ++i){
            _ffnn->setBeta(i, in[i]);
         }
      }
      
      void getVP(double *out){
         for (int i=0; i<_ffnn->getNBeta(); ++i){
            out[i] = _ffnn->getBeta(i);
         }
      }

      void samplingFunction(const double *x, double *out){
         /*
         Compute the sampling function proto value, used in getAcceptance()
         */
         _ffnn->setInput(1, x);
         _ffnn->FFPropagate();
         out[0] = pow(_ffnn->getOutput(1), 2);
      }

      double getAcceptance(){
         /*
         Compute the acceptance probability
         */
         if ((getProtoOld(0) == 0.) && (getProtoNew(0) != 0.)){
            return 1.;
         } else if ((getProtoOld(0) == 0.) && (getProtoNew(0) == 0.)) {
            return 0.;
         }
         
         return getProtoNew(0)/getProtoOld(0);
      }

      double d1(const int &i, const double *x){
         /*
         Compute:    d/dx_i log(Psi(x))
         */
         _ffnn->setInput(1, x);
         _ffnn->FFPropagate();
         return _ffnn->getFirstDerivative(1, i);
      }

      double d2(const int &i, const int &j, const double *x){
         /*
         Compute:    d^2/dx_i^2 log(Psi(x))
         */
         _ffnn->setInput(1, x);
         _ffnn->FFPropagate();
         return _ffnn->getSecondDerivative(1, i);
      }

      double vd1(const int &i, const double *x){
         /*
         Compute:    d/dalpha_i log(Psi(x))
         where alpha are the variational parameters, in our case an array of dimension 1: alpha = (b)
         */
         _ffnn->setInput(1, x);
         _ffnn->FFPropagate();
         return _ffnn->getVariationalFirstDerivative(1, i);
      }
};




int main(){
   using namespace std;

   // Declare some trial wave functions
   NeuralWaveFunction * psi = new NeuralWaveFunction(1, 30);
   
   // Store in a .txt file the values of the initial wf, so that it is possible to plot it
   cout << "Writing the plot file of the initial wave function in plot_init_wf.txt" << endl << endl;
   double * base_input = new double[psi->getFFNN()->getNInput()]; // no need to set it, since it is 1-dim
   const int input_i = 0;
   const int output_i = 1;
   const double min = -7.5;
   const double max = 7.5;
   const int npoints = 500;
   writePlotFile(psi->getFFNN(), base_input, input_i, output_i, min, max, npoints, "getOutput", "plot_init_wf.txt");
   
   
   // Declare an Hamiltonian
   // We use the harmonic oscillator with w=1 and w=2
   const double w1 = 1.;
   HarmonicOscillator1D1P * ham = new HarmonicOscillator1D1P(w1, psi);


   
   
   cout << endl << " - - - FFNN-WF FUNCTION OPTIMIZATION - - - " << endl << endl;
   
   VMC * vmc; // VMC object we will resuse
   const long E_NMC = 5000l; // MC samplings to use for computing the energy
   const long G_NMC = 10000l; // MC samplings to use for computing the energy gradient
   double * energy = new double[4]; // energy
   double * d_energy = new double[4]; // energy error bar
   double * vp = new double[psi->getNVP()];
      
   
   cout << "-> ham1:    w = " << w1 << endl << endl;
   vmc = new VMC(psi, ham);
   
   // set an integration range, because the NN might be completely delocalized
   double ** irange = new double*[1];
   irange[0] = new double[2];
   irange[0][0] = -7.5;
   irange[0][1] = 7.5;
   vmc->getMCI()->setIRange(irange);
   
   
   cout << "   Starting energy:" << endl;
   vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
   cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   
   cout << "   Optimization . . ." << endl;
   vmc->conjugateGradientOptimization(E_NMC, G_NMC);
   cout << "   . . . Done!" << endl << endl;
   
   cout << "   Optimized energy:" << endl;
   vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
   cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl << endl;
   
   
   // store in a .txt file the values of the optimised wf, so that it is possible to plot it
   cout << "Writing the plot file of the optimised wave function in plot_opt_wf.txt" << endl << endl;
   writePlotFile(psi->getFFNN(), base_input, input_i, output_i, min, max, npoints, "getOutput", "plot_opt_wf.txt");

   
   
   delete[] irange[0];
   delete[] irange;
   delete vmc;
   delete[] vp;
   delete[] d_energy;
   delete[] energy;
   delete ham;
   delete base_input;
   delete psi;

   

   return 0;
}


