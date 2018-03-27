#ifndef FUNCTION_WITH_VARIATIONAL_PARAMETERS
#define FUNCTION_WITH_VARIATIONAL_PARAMETERS



class FunctionWithVariationalParameters{
private:
    int _nvp;  //number of variational parameters involved

protected:
    void _setNVP(const int &nvp){_nvp = nvp;}

public:
    FunctionWithVariationalParameters(const int &nvp){
        _nvp = nvp;
    }
    virtual ~FunctionWithVariationalParameters(){
        _nvp = 0;
    }

    int getNVP(){return _nvp;}

    // --- methods that must be implemented
    virtual void getVP(double *vp) = 0;
    virtual void setVP(const double *vp) = 0;
};



#endif
