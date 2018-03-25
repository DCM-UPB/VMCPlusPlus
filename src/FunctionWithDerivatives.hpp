#ifndef FUNCTION_WITH_DERIVATIVES
#define FUNCTION_WITH_DERIVATIVES



class FunctionWithDerivatives{
protected:
    int _totalndim;
    int _nvp;

    bool _flag_d1;
    bool _flag_d2;
    bool _flag_vd1;
    bool _flag_d1vd1;
    bool _flag_d2vd1;

    double * _d1_divbywf;
    double * _d2_divbywf;
    double * _vd1_divbywf;
    double ** _d1vd1_divbywf;
    double ** _d2vd1_divbywf;



    // --- getters and setters for the derivatives (should be called only by child classes)
    // first derivative divided by the wf
    void _setD1DivByWF(const int &id1, const double &d1_divbywf){_d1_divbywf[id1] = d1_divbywf;}
    double * _getD1DivByWF(){return _d1_divbywf;}
    // second derivative divided by the wf
    void _setD2DivByWF(const int &id2, const double &d2_divbywf){_d2_divbywf[id2] = d2_divbywf;}
    double * _getD2DivByWF(){return _d2_divbywf;}
    // variational derivative divided by the wf
    void _setVD1DivByWF(const int &ivd1, const double &vd1_divbywf){_vd1_divbywf[ivd1] = vd1_divbywf;}
    double * _getVD1DivByWF(){return _vd1_divbywf;}
    // cross derivative: first derivative and first variational derivative divided by the wf
    void _setD1VD1DivByWF(const int &id1, const int &ivd1, const double &d1vd1_divbywf){_d1vd1_divbywf[id1][ivd1] = d1vd1_divbywf;}
    double ** _getD1VD1DivByWF(){return _d1vd1_divbywf;}
    // cross derivative: second derivative and first variational derivative divided by the wf
    void _setD2VD1DivByWF(const int &id2, const int &ivd1, const double &d2vd1_divbywf){_d2vd1_divbywf[id2][ivd1] = d2vd1_divbywf;}
    double ** _getD2VD1DivByWF(){return _d2vd1_divbywf;}


    void _allocateDerivativesMemory(const int &totalndim, const int &nvp);


public:
    FunctionWithDerivatives(const int &totalndim, const int &nvp, bool flag_d1, bool flag_d2, bool flag_vd1, bool flag_d1vd1, bool flag_d2vd1);
    ~FunctionWithDerivatives();


    // --- getters for the derivatives
    // first derivative divided by the wf
    bool hasD1(){return _flag_d1;}
    double getD1DivByWF(const int &id1){return _d1_divbywf[id1];}
    // second derivative divided by the wf
    bool hasD2(){return _flag_d2;}
    double getD2DivByWF(const int &id2){return _d2_divbywf[id2];}
    // variational derivative divided by the wf
    bool hasVD1(){return _flag_vd1;}
    double getVD1DivByWF(const int &ivd1){return _vd1_divbywf[ivd1];}
    // cross derivative: first derivative and first variational derivative divided by the wf
    bool hasD1VD1(){return _flag_d1vd1;}
    double getD1VD1DivByWF(const int &id1, const int &ivd1){return _d1vd1_divbywf[id1][ivd1];}
    // cross derivative: second derivative and first variational derivative divided by the wf
    bool hasD2VD1(){return _flag_d2vd1;}
    double getD2VD1DivByWF(const int &id2, const int &ivd1){return _d2vd1_divbywf[id2][ivd1];}



    // --- computation of the derivatives
    // When called, this method computes all the internal values, such as the derivatives,
    // and stored them internally, ready to be accessed with the getters methods.
    // It requires the positions as input
    virtual void computeAllDerivatives(const double *x) = 0;    // --- MUST BE IMPLEMENTED


};



#endif
