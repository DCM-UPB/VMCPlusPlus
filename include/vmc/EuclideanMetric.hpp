#ifndef VMC_EUCLIDEANMETRIC_HPP
#define VMC_EUCLIDEANMETRIC_HPP

#include "vmc/Metric.hpp"

namespace vmc
{

class EuclideanMetric: public Metric
{
public:
    explicit EuclideanMetric(const int &nspacedim): Metric(nspacedim){}
    ~EuclideanMetric() override= default;

    double dist(const double * r1, const double * r2) override;

    void distD1(const double * r1, const double * r2, double * out) override;

    void distD2(const double * r1, const double * r2, double * out) override;
};
} // namespace vmc

#endif
