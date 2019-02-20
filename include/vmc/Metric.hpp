#ifndef METRIC
#define METRIC


class Metric
{
private:
    int _nspacedim;

public:
    Metric(const int &nspacedim){
        _nspacedim = nspacedim;
    }
    virtual ~Metric(){}


    int getNSpaceDim(){return _nspacedim;}


    // --- Methods that must be implemented
    virtual double dist(const double * r1, const double * r2) = 0;

    virtual void distD1(const double * r1, const double * r2, double * out) = 0;

    virtual void distD2(const double * r1, const double * r2, double * out) = 0;
};


#endif
