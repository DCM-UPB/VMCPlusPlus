#ifndef FUNCTION_WITH_VARIATIONAL_PARAMETERS
#define FUNCTION_WITH_VARIATIONAL_PARAMETERS



class FunctionWithVariationalParameters{
private:
    int _nvp;  //number of variational parameters involved
    double * _vp;

public:
    FunctionWithVariationalParameters(const int &nvp);
    ~FunctionWithVariationalParameters();

    int getNVP();
    void setNVP(const int &nvp);

    void getVP(double *vp);
    double getVP(const int &i);

    void setVP(const double *vp);
    void setVP(const int &i, const double &vp);

    virtual void actAfterVPChange(const int &i, const double &vp) = 0;
};



#endif
