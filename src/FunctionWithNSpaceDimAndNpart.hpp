#ifndef FUNCTION_WITH_NSPACEDIM_AND_NPART
#define FUNCTION_WITH_NSPACEDIM_AND_NPART



class FunctionWithNSpaceDimAndNPart{
private:
    int _nspacedim;
    int _npart;

public:
    FunctionWithNSpaceDimAndNPart(const int &nspacedim, const int &npart);
    ~FunctionWithNSpaceDimAndNPart();

    int getNSpaceDim();
    int getTotalNDim();
    int getNPart();
};



#endif
