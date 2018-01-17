#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>

#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"


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
Gaussian Trial Wave Function for 1 particle in 1 dimension, that uses two variational parameters: x0 and v (variance).
Psi  =  exp( -0.5*(x-x0)^2/v )
Notice that the corresponding probability density (sampling function) is Psi^2.
*/
class Gaussian1D1POrbital: public WaveFunction{
   protected:
  double _x0, _v, _niv;

   public:
      Gaussian1D1POrbital(const double x0, const double v): 
        WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, 2 /*num variational parameters*/) {_x0=x0; _v=v; _niv = -1./v;}

      void setVP(const double *in){
         _x0=in[0];
         _v=in[1];
         _niv = -1./in[1];
      }

      void getVP(double *out){
         out[0]=_x0; out[1]=_v;
      }

      void samplingFunction(const double *x, double *out){
         /*
         Compute the sampling function proto value, used in getAcceptance()
         */
        *out = _niv*(x[0]-_x0)*(x[0]-_x0);
      }

      double getAcceptance(){
         /*
         Compute the acceptance probability
         */
         return exp(getProtoNew(0)-getProtoOld(0));
      }

      double d1(const int &i, const double *x){
         /*
         Compute:    d/dx_i log(Psi(x))
         */
         return (_niv*(x[0]-_x0) );
      }

      double d2(const int &i, const int &j, const double *x){
         /*
         Compute:    d^2/dx_i^2 log(Psi(x))
         */
         return ( _niv + (_niv*(x[0]-_x0))*(_niv*(x[0]-_x0)) ) ;
      }

      double vd1(const int &i, const double *x){
         /*
         Compute:    d/dalpha_i log(Psi(x))
         where alpha are the variational parameters, in our case an array of dimension 1: alpha = (b)
         */
         if (i==0){
            return (-_niv*(x[0]-_x0));
         } else if (i==1){
            return (-(x[0]-_x0)*(x[0]-_x0));
         } else{
            using namespace std;
            cout << "ERROR vd1 Gaussian1D1POrbital! " << endl;
            return 0.;
         }
      }
};


class FFNNGaussian1D1P: public WaveFunction{
protected:
  Gaussian1D1POrbital * _gauss;
  FeedForwardNeuralNetwork * _ffnn;

public:
  FFNNGaussian1D1P(const double x0, const double v, const int nh1, const int nh2):
    WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, nh1-1 + ((nh2>1)? nh2-1 : 0) /*num variational parameters*/) {
    using namespace std;
    _gauss = new Gaussian1D1POrbital(x0, v);
    _ffnn = new FeedForwardNeuralNetwork(2, nh1, 2);
    if (nh2>1) {_ffnn->pushHiddenLayer(nh2);}

    cout << "Creating wavefunction based on a Feed Forward Artificial Neural Network (FFANN)" << endl;
    cin.ignore();

    cout << "Graphically it looks like this" << endl;
    cin.ignore();
    printFFNNStructure(_ffnn);
    cout << endl << endl;
    cin.ignore();

    cout << "Connecting the FFNN..." << endl;
    cin.ignore();

    // NON I/O CODE
    _ffnn->connectFFNN();
    //

    cout << endl << endl;


    cout << "Adding derivatives substrates..." << endl;
    cin.ignore();

    // NON I/O CODE
    _ffnn->addFirstDerivativeSubstrate();
    _ffnn->addSecondDerivativeSubstrate();
    //

    cout << endl << endl;


    cout << "Setting the input to 0..." << endl;
    cin.ignore();

    int ninput = 1;
    double * input = new double[ninput];
    input[0] = 0;

    // NON I/O CODE
    _ffnn->setInput(ninput, input);
    //
    cout << "Done! Now the NN values look like this:";
    cin.ignore();

    printFFNNValues(_ffnn);

    cin.ignore();
    cout << endl << endl;

    cout << "Propagating..." << endl;
    cin.ignore();

    // NON I/O CODE
    _ffnn->FFPropagate();
    //

    cout << "Done! Now the NN values look like this:";
    cin.ignore();

    printFFNNValues(_ffnn);
    cout << endl;
    cin.ignore();

    cout << "The output values are ";
    cout << _ffnn->getOutput(1) << endl;
    cin.ignore();

    cout << "The first derivative with respect to the input value is:";
    cin.ignore();
    cout << "1st output (unit 2 of the output layer): " << _ffnn->getFirstDerivative(1, 0) << endl;
    cin.ignore();

    cout << "The second derivative with respect to the input value is:";
    cin.ignore();
    cout << "1st output (unit 2 of the output layer): " << _ffnn->getSecondDerivative(1, 0);
    cin.ignore();

    cout << endl << endl;

  }

  void setX(const double x){
    using namespace std;
    double *xarr = new double(1);
    xarr[0] = x;

    //cout << "Setting the input to x..." << endl;
    //cin.ignore();

    cout << "The input we want to set is: " << xarr[0];
    //cin.ignore();

    // NON I/O CODE
    _ffnn->setInput(1, xarr);
    //
    //cout << "Done! Now the NN values look like this:";
    //cin.ignore();

    //printFFNNValues(_ffnn);

    //cin.ignore();
    //cout << endl << endl;

    //cout << "Propagating..." << endl;
    //cin.ignore();

    // NON I/O CODE
    _ffnn->FFPropagate();
    //

    // cout << "Done! Now the NN values look like this:";
    //cin.ignore();

    //printFFNNValues(_ffnn);
    //cout << endl;
    //cin.ignore();
  }

  double Psi(const double x){
    using namespace std;
    double out;

    setX(x);

    out = _ffnn->getOutput(1);

    cout << "The output value is " << out << endl;
    //cin.ignore();

    return out;
  }

  void setVP(const double *in){
  }

  void getVP(double *out){
  }

  void samplingFunction(const double *x, double *out){
    /*
      Compute the sampling function proto value, used in getAcceptance()
    */
    double psix = Psi(x[0]);
    out[0] = log(psix*psix);
  }

  double getAcceptance(){
    /*
      Compute the acceptance probability
    */
    return exp(getProtoNew(0)-getProtoOld(0));
  }

  double d1(const int &i, const double *x){
    using namespace std;
    /*
      Compute:    d/dx_i log(Psi(x))
    */

    double psix = Psi(x[0]);
    double out = _ffnn->getFirstDerivative(1, 0) / psix;

    //cout << "The first derivative with respect to the input value is:";
    //cin.ignore();
    //cout << "1st output (unit 2 of the output layer): " << out << endl;
    //cin.ignore();

    return out;
  }

  double d2(const int &i, const int &j, const double *x){
    using namespace std;
    /*
      Compute:    d^2/dx_i^2 log(Psi(x))
    */

    double psix = Psi(x[0]);
    double psix2 = psix*psix;
    double psid1 = _ffnn->getFirstDerivative(1, 0);
    double psid2 = _ffnn->getSecondDerivative(1, 0);
    double out = (psid2*psix - psid1*psid1) / psix2;

    //cout << "The second derivative with respect to the input value is:";
    //cin.ignore();
    //cout << "1st output (unit 2 of the output layer): " << out << endl;;
    //cin.ignore();


    return out;
  }

  double vd1(const int &i, const double *x){
    /*
      Compute:    d/dalpha_i log(Psi(x))
      where alpha are the variational parameters, in our case an array of dimension 1: alpha = (b)
    */
    if (i==0){
      return 0;
    } else if (i==1){
      return 0;
    } else{
      using namespace std;
      cout << "ERROR vd1 FFNNGaussian1D1P! " << endl;
      return 0.;
    }
  }
};


int main() {
    using namespace std;

    int nl, nh1, nh2;

    cout << "Let's start by creating a Feed Forward Artificial Neural Netowrk (FFANN)" << endl;
    cout << "========================================================================" << endl;
    cin.ignore();

    cout << "How many units should the first hidden layer(s) have? ";
    cin >> nh1;
    cout << "How many units should the first hidden layer(s) have? (<=1 for no second hidden layer)";
    cin >> nh2;

    nl = (nh2>1)? 4:3; 
    cout << "We generate a FFANN with " << nl << " layers and 2, " << nh1;
    if (nh2>1) { cout << ", " << nh2;}
    cout << ", 2 units respectively" << endl;
    cout << "========================================================================" << endl << endl;
    cin.ignore();


    // NON I/O CODE
    FFNNGaussian1D1P * nnwf = new FFNNGaussian1D1P(0, 1, nh1, nh2);
    //
    cout << "Psi(1):" << endl;
    cout << nnwf->Psi(1) << endl;;
    cout << "acceptance: " <<  nnwf->getAcceptance() << endl;
    cin.ignore();


   // Declare some trial wave functions
    Gaussian1D1POrbital * psi1 = new Gaussian1D1POrbital(0., 1);
    Gaussian1D1POrbital * psi2 = new Gaussian1D1POrbital(0., 2);
    Gaussian1D1POrbital * psi3 = new Gaussian1D1POrbital(1., 1);
    Gaussian1D1POrbital * psi4 = new Gaussian1D1POrbital(1., 2.);

   // Declare an Hamiltonian for each wave function (keep in mind that the kinetic energy is strictly bound to it)
   // We use the harmonic oscillator with w=1
   HarmonicOscillator1D1P * ham1 = new HarmonicOscillator1D1P(1., psi1);
   HarmonicOscillator1D1P * ham2 = new HarmonicOscillator1D1P(1., psi2);
   HarmonicOscillator1D1P * ham3 = new HarmonicOscillator1D1P(1., psi3);
   HarmonicOscillator1D1P * ham4 = new HarmonicOscillator1D1P(1., psi4);

   HarmonicOscillator1D1P * hamnn = new HarmonicOscillator1D1P(1., nnwf);


   cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;

   VMC * vmc; // VMC object we will resuse
   const long E_NMC = 100000l; // MC samplings to use for computing the energy
   double * energy = new double[4]; // energy
   double * d_energy = new double[4]; // energy error bar

   // Case 1
   cout << "-> psi1: " << endl;
   vmc = new VMC(psi1, ham1, 0, 0); // ENMC and GNMC do not need to be set since we don't optimize the wave function
   vmc->computeEnergy(E_NMC, energy, d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

   // Case 2
   cout << "-> psi2: " << endl;
   delete vmc;
   vmc = new VMC(psi2, ham2, 0, 0);
   vmc->computeEnergy(E_NMC, energy, d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

   // Case 3
   cout << "-> psi3: " << endl;
   delete vmc;
   vmc = new VMC(psi3, ham3, 0, 0);
   vmc->computeEnergy(E_NMC, energy, d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

   // Case 4
   cout << "-> psi4: " << endl;
   delete vmc;
   vmc = new VMC(psi4, ham4, 0, 0);
   vmc->computeEnergy(E_NMC, energy, d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

   // Case NN
   cout << "-> psinn: " << endl;
   delete vmc;
   vmc = new VMC(nnwf, hamnn, 0, 0);
   vmc->computeEnergy(2l, energy, d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;


   //deallocate everything
   delete vmc;

   delete psi1;
   delete psi2;
   delete psi3;
   delete psi4;

   delete ham1;
   delete ham2;
   delete ham3;
   delete ham4;
   delete hamnn;

   delete[] d_energy;
   delete[] energy;

    // end
    return 0;
}
