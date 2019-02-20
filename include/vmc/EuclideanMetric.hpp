#ifndef EUCLIDEAN_METRIC
#define EUCLIDEAN_METRIC

#include "vmc/Metric.hpp"


class EuclideanMetric: public Metric{

public:
    explicit EuclideanMetric(const int &nspacedim): Metric(nspacedim){}
    ~EuclideanMetric(){}

    double dist(const double * r1, const double * r2);

    void distD1(const double * r1, const double * r2, double * out);

    void distD2(const double * r1, const double * r2, double * out);

};


#endif
