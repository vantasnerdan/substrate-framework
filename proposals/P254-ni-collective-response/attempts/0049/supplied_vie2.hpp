// Supplied contour-Hermitian R/< boundary; pinned upstream is the integrator.
#ifndef P254_SUPPLIED_VIE2_HPP
#define P254_SUPPLIED_VIE2_HPP
#include "../0040/nessi_gaussian_ivp.hpp"

namespace p254::supplied {
inline double uniform_step(const std::vector<double>& times, int order) {
    if (order<1 || order>5 || times.size()<=size_t(order))
        throw std::invalid_argument("complete startup block 0..order required; order 1..5");
    const double dt=times[1]-times[0];
    if (times[0]!=0. || !(dt>0.) || !std::isfinite(dt))
        throw std::invalid_argument("uniform real time must start at zero with positive step");
    for (size_t n=0;n<times.size();++n)
        if (!std::isfinite(times[n]) || std::abs(times[n]-n*dt)>1.e-12*std::max(1.,std::abs(n*dt)))
            throw std::invalid_argument("nonuniform sample times unsupported; no interpolation");
    return dt;
}

// All objects have the same literal statistics and zero initialized imaginary
// branch. Source/solution must be representable as one contour-Hermitian object.
// These functions never infer arbitrary independent reverse functional legs.
class Solver {
    int order_;
    double dt_, beta_;
    Green f_, fcc_;
public:
    Solver(const std::vector<double>& times, int dimension, int statistics,
           int order, double beta=1.)
        : order_(order), dt_(uniform_step(times,order)), beta_(beta),
          f_(times.size()-1,std::max(2*order+2,12),dimension,statistics),
          fcc_(times.size()-1,std::max(2*order+2,12),dimension,statistics) {
        if (dimension<1 || (statistics!=-1 && statistics!=1) || !(beta>0.) || !std::isfinite(beta))
            throw std::invalid_argument("invalid dimension/statistics/inert beta");
        clear(f_); clear(fcc_);
    }
    void check(const Green& g) const {
        if (g.nt()!=f_.nt() || g.ntau()!=f_.ntau() || g.size1()!=f_.size1() || g.sig()!=f_.sig())
            throw std::invalid_argument("source/kernel/solution layout or statistics mismatch");
    }
    // Refresh after the caller imports its current common-functional kernels.
    // Early slices read the COMPLETE 0..order source/kernel block.
    void refresh(int n, Green& source, Green& kernel) {
        check(source); check(kernel);
        if (n<0 || n>f_.nt()) throw std::out_of_range("refresh time slice");
        auto& integrator=integration::I<double>(order_);
        cntr::convolution_timestep(n,f_,source,kernel,integrator,beta_,dt_);
        cntr::convolution_timestep(n,fcc_,kernel,source,integrator,beta_,dt_);
        f_.smul(n,-1.); fcc_.smul(n,-1.);
    }
    void start(Green& result, Green& source, Green& kernel) {
        check(result);
        for (int n=0;n<=order_;++n) refresh(n,source,kernel);
        cntr::vie2_start(result,f_,fcc_,source,integration::I<double>(order_),beta_,dt_);
    }
    void step(int n, Green& result, Green& source, Green& kernel) {
        check(result);
        if (n<=order_ || n>f_.nt()) throw std::out_of_range("step must follow complete startup");
        refresh(n,source,kernel);
        cntr::vie2_timestep(n,result,f_,fcc_,source,integration::I<double>(order_),beta_,dt_);
    }
    const Green& forward_kernel() const { return f_; }
    const Green& reverse_kernel() const { return fcc_; }
};
}
#endif
