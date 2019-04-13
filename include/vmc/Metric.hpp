#ifndef VMC_METRIC_HPP
#define VMC_METRIC_HPP

namespace vmc
{

class Metric
{
private:
    const int _nspacedim;

public:
    explicit Metric(int nspacedim): _nspacedim(nspacedim) {}

    virtual ~Metric() = default;


    int getNSpaceDim() { return _nspacedim; }


    // --- Methods that must be implemented
    virtual double dist(const double * r1, const double * r2) = 0;

    virtual void distD1(const double * r1, const double * r2, double * out) = 0;

    virtual void distD2(const double * r1, const double * r2, double * out) = 0;
};
} // namespace vmc

#endif
