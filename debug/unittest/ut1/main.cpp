#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "FFNNWaveFunction.hpp"
#include "FeedForwardNeuralNetwork.hpp"
#include "GaussianActivationFunction.hpp"
#include "IdentityActivationFunction.hpp"
#include "VMC.hpp"
#include "PrintUtilities.hpp"

#include <iostream>
#include <assert.h>
#include <math.h>
#include <stdexcept>


/*

In this unit test we test that the FFNNWaveFunction works as expected.
In particular we want to be sure that all the derivatives are computed properly.

To accomplish this, we are following this approach.
We build a NN with the following structure:

    id_            id_            id_

    id_     p1     gss     +0.00  id_
            p2             +1.00

where p1 nad p2 are two variational parameters.
This means that the NN is the following function:

    nn(x) = exp( - ( p1 + p2*x )^2 )

i.e. a gaussian. We write a WaveFunction that has the same structure:

    gauss(x) = exp( - b ( x - a )^2 )

There is therefore the following mapping:

    p1  = - sqrt(b) * a
    p2 = sqrt(b)

Setting the parameters respecting the mapping, the two functions and therefore should give the same result.
We then checked that the FFNNWaveFunction works as expected with two tests:

1. Computing the variational energy with the two wave function must return the same result

2. Computing the variational derivatives must result into the same results, asides for some costants.
   In particular:

      d(nn)/da = d(p1)/da * d(nn)/dp1 = -sqrt(b) * d(nn)/dp1

      d(nn)/db = d(p1)/db * d(nn)/dp1 + d(p2)/db * d(nn)/dp2 = ( -a / (2 * sqrt(b)) ) * d(nn)/dp1 + ( 1 / (2 * sqrt(b)) ) * d(nn)/dp2


*/





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
Trial Wave Function for 1 particle in 1 dimension, that uses two variational parameters: a and b.
   Psi  =  exp( -b * (x-a)^2 )
Notice that the corresponding probability density (sampling function) is Psi^2.
*/
class QuadrExponential1D1POrbital: public WaveFunction{
   protected:
      double _a, _b;

   public:
      QuadrExponential1D1POrbital(const double a, const double b):
         WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, 2 /*num variational parameters*/) {_a=a; _b=b;}

      void setVP(const double *in){
         _a=in[0];
         _b=in[1];
      }

      void getVP(double *out){
         out[0]=_a; out[1]=_b;
      }

      void samplingFunction(const double *x, double *out){
         /*
         Compute the sampling function proto value, used in getAcceptance()
         */
         *out = -2.*(_b*(x[0]-_a)*(x[0]-_a));
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
         return (-2.*_b*(x[0]-_a) );
      }

      double d2(const int &i, const double *x){
         /*
         Compute:    d^2/dx_i^2 log(Psi(x))
         */
         return ( -2.*_b + (-2.*_b*(x[0]-_a))*(-2.*_b*(x[0]-_a)) ) ;
      }

      double vd1(const int &i, const double *x){
         /*
         Compute:    d/dalpha_i log(Psi(x))
         where alpha are the variational parameters, in our case an array of dimension 1: alpha = (b)
         */
         if (i==0){
            return (2.*_b*(x[0]-_a));
         } else if (i==1){
            return (-(x[0]-_a)*(x[0]-_a));
         } else{
            throw std::range_error( " the index i for QuadrExponential1D1POrbital.vd1() can be only 0 or 1" );
         }
      }
};



int main(){
    using namespace std;

    // parameters
    const long Nmc = 100000l;
    const double TINY = 0.000001;

    // variational parameters
    double a = 0.37;
    double sqrtb = 1.18;
    double b = sqrtb * sqrtb;
    double p1 = - sqrtb * a;
    double p2 = sqrtb;




    // create the NN with the right structure and parameters
    FeedForwardNeuralNetwork * ffnn = new FeedForwardNeuralNetwork(2, 2, 2);
    GaussianActivationFunction * gss_actf = new GaussianActivationFunction();
    IdentityActivationFunction * id_actf = new IdentityActivationFunction();
    ffnn->setLayerActivationFunction(1, gss_actf);
    ffnn->setLayerActivationFunction(2, id_actf);
    ffnn->connectFFNN();
    ffnn->addFirstDerivativeSubstrate();
    ffnn->addSecondDerivativeSubstrate();
    ffnn->addVariationalFirstDerivativeSubstrate();

    ffnn->setBeta(0, p1);
    ffnn->setBeta(1, p2);
    ffnn->setBeta(2, 0.);
    ffnn->setBeta(3, 1.);

    // printFFNNStructure(ffnn);
    // printFFNNStructureWithBeta(ffnn);


    // NN wave function
    FFNNWaveFunction * psi = new FFNNWaveFunction(1, 1, ffnn);

    // gaussian wave function
    QuadrExponential1D1POrbital * phi = new QuadrExponential1D1POrbital(a, b);

    // Hamiltonian
    HarmonicOscillator1D1P * ham = new HarmonicOscillator1D1P(1., psi);






    // --- Check that the energies are the same
    VMC * vmc = new VMC(psi, ham);

    double energy[4];
    double d_energy[4];
    vmc->computeVariationalEnergy(Nmc, energy, d_energy);
    // cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    // cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    // cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    // cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;


    double energy_check[4];
    double d_energy_check[4];
    VMC * vmc_check = new VMC(phi, ham);
    vmc_check->computeVariationalEnergy(Nmc, energy_check, d_energy_check);
    // cout << "       Total Energy        = " << energy_check[0] << " +- " << d_energy_check[0] << endl;
    // cout << "       Potential Energy    = " << energy_check[1] << " +- " << d_energy_check[1] << endl;
    // cout << "       Kinetic (PB) Energy = " << energy_check[2] << " +- " << d_energy_check[2] << endl;
    // cout << "       Kinetic (JF) Energy = " << energy_check[3] << " +- " << d_energy_check[3] << endl << endl;

    for (int i=0; i<4; ++i){
        assert(abs(energy[i]-energy_check[i]) < 2.*(d_energy[i]+d_energy_check[i]));
    }






    // --- Check the variational derivatives
    const double dx=0.2;
    double x = -1.;
    for (int i=0; i<10; ++i){
        x = x + dx;

        const double dda_phi = phi->vd1(0, &x);
        const double ddb_phi = phi->vd1(1, &x);

        const double ddp1_psi = psi->vd1(0, &x);
        const double ddp2_psi = psi->vd1(1, &x);

        // cout << dda_phi << " == " << - ddp1_psi * sqrtb << " ? " << endl;
        assert( abs(dda_phi - ( - ddp1_psi * sqrtb )) < TINY );

        // cout << ddb_phi << " == " << ddp2_psi / (2. * sqrtb) - ddp1_psi * a / (2. * sqrtb) << " ? " << endl;
        assert( abs(ddb_phi - ( ddp2_psi / (2. * sqrtb) - ddp1_psi * a / (2. * sqrtb) ) ) < TINY );
    }





    // free resources
    delete phi;
    delete vmc;
    delete ham;
    delete psi;
    delete ffnn;


    return 0;
}
