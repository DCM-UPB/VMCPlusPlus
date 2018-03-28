#ifndef FUNCTION_WITH_NSPACEDIM_AND_NPART
#define FUNCTION_WITH_NSPACEDIM_AND_NPART

#include <stdexcept>



class FunctionWithNSpaceDimAndNPart{
private:
    int _nspacedim;
    int _npart;

public:
    FunctionWithNSpaceDimAndNPart(const int &nspacedim, const int &npart){
        if (nspacedim <= 0){
            throw std::invalid_argument( "nspacedim must be strictly greater than zero" );
        }
        if (npart <= 0){
            throw std::invalid_argument( "npart must be strictly greater than zero" );
        }
        
        _nspacedim = nspacedim;
        _npart = npart;
    }
    ~FunctionWithNSpaceDimAndNPart(){
        _nspacedim = 0;
        _npart = 0;
    }

    int getNSpaceDim(){return _nspacedim;}
    int getTotalNDim(){return _nspacedim * _npart;}
    int getNPart(){return _npart;}
};



#endif
