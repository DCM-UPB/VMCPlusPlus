#ifndef VMC_WAVEFUNCTION_HPP
#define VMC_WAVEFUNCTION_HPP

#include "mci/CallBackOnMoveInterface.hpp"
#include "mci/SamplingFunctionInterface.hpp"

#include <iostream>

namespace vmc
{
/*
IMPLEMENTATIONS OF THIS INTERFACE MUST INCLUDE:

    - void setVP(const double *vp)
            set the variational parameters

    - void getVP(double *vp)
            get the variational parameters

    - void protoFunction(const double * in, double * out)
            heritage from mci::SamplingFunctionInterface, uses Psi^2

    - double acceptanceFunction(const double * protoold, const double * protonew)
            heritage from mci::SamplingFunctionInterface

    - double computeWFValue(const double * protov)
            compute the true wave function value from given proto values

    - void computeAllDerivatives(const double *x)
            use the setters for derivatives values (setD1DivByWF, setD2DivByWF, etc.)

*/
class WaveFunction: public mci::SamplingFunctionInterface
{
protected:
    const int _nspacedim;
    const int _npart;
    int _nvp;  //number of variational parameters involved

    const bool _flag_vd1;
    const bool _flag_d1vd1;
    const bool _flag_d2vd1;

    double * _d1_divbywf;
    double * _d2_divbywf;
    double * _vd1_divbywf;
    double ** _d1vd1_divbywf;
    double ** _d2vd1_divbywf;


    void _allocateVariationalDerivativesMemory();

    // --- getters and setters for the derivatives
    // first derivative divided by the wf
    void _setD1DivByWF(int id1, double d1_divbywf) { _d1_divbywf[id1] = d1_divbywf; }
    double * _getD1DivByWF() const { return _d1_divbywf; }
    // second derivative divided by the wf
    void _setD2DivByWF(int id2, double d2_divbywf) { _d2_divbywf[id2] = d2_divbywf; }
    double * _getD2DivByWF() const { return _d2_divbywf; }
    // variational derivative divided by the wf
    void _setVD1DivByWF(int ivd1, double vd1_divbywf) { _vd1_divbywf[ivd1] = vd1_divbywf; }
    double * _getVD1DivByWF() const { return _vd1_divbywf; }
    // cross derivative: first derivative and first variational derivative divided by the wf
    void _setD1VD1DivByWF(int id1, int ivd1, double d1vd1_divbywf) { _d1vd1_divbywf[id1][ivd1] = d1vd1_divbywf; }
    double ** _getD1VD1DivByWF() const { return _d1vd1_divbywf; }
    // cross derivative: second derivative and first variational derivative divided by the wf
    void _setD2VD1DivByWF(int id2, int ivd1, double d2vd1_divbywf) { _d2vd1_divbywf[id2][ivd1] = d2vd1_divbywf; }
    double ** _getD2VD1DivByWF() const { return _d2vd1_divbywf; }

public:
    WaveFunction(int nspacedim, int npart, int ncomp/*defines number of proto values*/,
                 int nvp, bool flag_vd1 = false, bool flag_d1vd1 = false, bool flag_d2vd1 = false);
    ~WaveFunction() override;

    int getNSpaceDim() const { return _nspacedim; }
    int getTotalNDim() const { return this->getNDim(); }
    int getNPart() const { return _npart; }
    int getNVP() const { return _nvp; }


    // --- interface for manipulating the variational parameters
    void setNVP(int nvp);
    virtual void setVP(const double * vp) = 0;    // --- MUST BE IMPLEMENTED
    virtual void getVP(double * vp) const = 0;    // --- MUST BE IMPLEMENTED


    // --- computation of the derivatives
    // When called, this method computes all the internal values, such as the derivatives,
    // and stored them internally, ready to be accessed with the getters methods.
    // It requires the positions as input
    virtual void computeAllDerivatives(const double * x) = 0;    // --- MUST BE IMPLEMENTED

    // --- computation of the wavefunction value
    // As the sampling function routine doesn't provide the actual
    // wavefunction value, you have to provide a method to reconstruct
    // the wactual wavefunction value (not squared) from the sampling
    // function values, in case it is required e.g. by a wrapper.
    // This method is also used to provide MCI's samplingFunction method.
    virtual double computeWFValue(const double * protovalues) const = 0;    // --- MUST BE IMPLEMENTED

    double samplingFunction(const double protovalues[]) const final
    { // mainly for use by certain MCI trial moves
        const double wfval = this->computeWFValue(protovalues);
        return wfval*wfval; // the sampling function is Psi^2
    }

    // --- getters and setters for the derivatives
    // first derivative divided by the wf
    double getD1DivByWF(int id1) const { return _d1_divbywf[id1]; }
    // second derivative divided by the wf
    double getD2DivByWF(int id2) const { return _d2_divbywf[id2]; }
    // variational derivative divided by the wf
    bool hasVD1() const { return _flag_vd1; }
    double getVD1DivByWF(int ivd1) const { return _vd1_divbywf[ivd1]; }
    // cross derivative: first derivative and first variational derivative divided by the wf
    bool hasD1VD1() const { return _flag_d1vd1; }
    double getD1VD1DivByWF(int id1, int ivd1) const { return _d1vd1_divbywf[id1][ivd1]; }
    // cross derivative: second derivative and first variational derivative divided by the wf
    bool hasD2VD1() const { return _flag_d2vd1; }
    double getD2VD1DivByWF(int id2, int ivd1) const { return _d2vd1_divbywf[id2][ivd1]; }
};
} // namespace vmc

#endif
