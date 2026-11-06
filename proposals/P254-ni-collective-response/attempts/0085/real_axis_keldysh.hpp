// Result-only implementation of A54's finite-step R/< operator.
#ifndef P254_A85_REAL_AXIS_KELDYSH_HPP
#define P254_A85_REAL_AXIS_KELDYSH_HPP
#include "../0049/component_io.hpp"
#include <limits>

namespace p254::real_axis {
// The P254VIE1 representation contains ONLY R and <. There is no imaginary
// input, equilibrium preparation, or mixed-channel fallback in this type.
// ntau is inert metadata required by the upstream REAL-AXIS helper interface;
// no Matsubara/mixed storage or accessors exist. Missing real channels are an
// error, not implicit zeros. Callers select the nonzero Langreth terms explicitly.
class History {
    int nt_, dimension_, statistics_, ntau_;
    size_t entries_;
    std::vector<Complex> ret_, les_;
    size_t offset(int n, int m) const {
        if (m<0 || n<m || n>nt_) throw std::out_of_range("real-axis triangle index");
        return (size_t(n)*(size_t(n)+1)/2+size_t(m))*size_t(element_size());
    }
public:
    History(int nt, int dimension, int statistics, int order, bool ret, bool les)
        : nt_(nt), dimension_(dimension), statistics_(statistics),
          ntau_(std::max(2*order+2,12)), entries_(0) {
        if (order<1 || order>5 || nt<order || dimension<1 ||
            (statistics!=-1 && statistics!=1) ||
            dimension>std::numeric_limits<int>::max()/dimension)
            throw std::invalid_argument("invalid real-axis history layout");
        const size_t nodes=size_t(nt)+1, cells=nodes*(nodes+1)/2;
        if (cells>std::vector<Complex>().max_size()/size_t(element_size()))
            throw std::invalid_argument("real-axis history exceeds native allocation representation");
        entries_=cells*size_t(element_size());
        if (ret) ret_.resize(entries_);
        if (les) les_.resize(entries_);
    }
    History(const History&)=delete;
    History& operator=(const History&)=delete;
    int nt() const { return nt_; }
    int ntau() const { return ntau_; }
    int size1() const { return dimension_; }
    int element_size() const { return dimension_*dimension_; }
    int sig() const { return statistics_; }
    Complex* retptr(int n,int m) {
        if (ret_.empty()) throw std::logic_error("absent retarded channel");
        return ret_.data()+offset(n,m);
    }
    Complex* lesptr(int m,int n) {
        if (les_.empty()) throw std::logic_error("absent lesser channel");
        return les_.data()+offset(n,m);
    }
    void negate_ret() { for (auto& x:ret_) x=-x; }
};

inline void read_component(std::istream& in, History& g, bool lesser) {
    Matrix diagonal(g.size1(),g.size1());
    for (int n=0;n<=g.nt();++n) for (int m=0;m<=n;++m) {
        auto* p=lesser ? g.lesptr(m,n) : g.retptr(n,m);
        for (int r=0;r<g.size1();++r) for (int s=0;s<g.size1();++s) {
            double re,im; supplied::read_value(in,re); supplied::read_value(in,im);
            if (!std::isfinite(re) || !std::isfinite(im))
                throw std::runtime_error("nonfinite supplied component");
            p[r*g.size1()+s]={re,im};
            if (lesser && n==m) diagonal(r,s)={re,im};
        }
        if (lesser && n==m &&
            (diagonal+diagonal.adjoint()).norm()>1.e-10*(1.+diagonal.norm()))
            throw std::runtime_error("lesser diagonal not contour-Hermitian; no repair");
    }
}
inline void write_component(std::ostream& out, History& g, bool lesser) {
    for (int n=0;n<=g.nt();++n) for (int m=0;m<=n;++m) {
        const auto* p=lesser ? g.lesptr(m,n) : g.retptr(n,m);
        for (int z=0;z<g.element_size();++z) {
            if (!std::isfinite(p[z].real()) || !std::isfinite(p[z].imag()))
                throw std::runtime_error("nonfinite solver output");
            supplied::write_value(out,p[z].real()); supplied::write_value(out,p[z].imag());
        }
    }
}

class Solver {
    int order_, dimension_, statistics_, nt_;
    double dt_, beta_;
    integration::Integrator<double>& integrator_;
    // One reusable timestep, including the complete startup extension 0..order.
    std::vector<Complex> work_;
    void check(History& g) const {
        if (g.nt()!=nt_ || g.size1()!=dimension_ || g.sig()!=statistics_ ||
            g.ntau()!=std::max(2*order_+2,12))
            throw std::invalid_argument("source/kernel/solution layout or statistics mismatch");
    }
    void ret_product(History& out,History& a,History& acc,History& b,History& bcc) {
        for (int n=0;n<=nt_;++n)
            cntr::convolution_timestep_ret<double,History,LARGESIZE>(
                n,out,a,acc,b,bcc,integrator_,dt_);
    }
    // Exactly one nonzero term of the old full-contour product. Upstream
    // helpers retain their original polynomial/Gregory startup continuations,
    // including distinct companions for noncommuting products.
    void lesser_slice(int n, History& out,History& a,History& acc,
                      History& b,History& bcc,bool retles) {
        const size_t count=size_t(std::max(n,order_)+1)*out.element_size();
        std::fill(work_.begin(),work_.begin()+count,Complex(0.));
        if (retles)
            cntr::convolution_timestep_les_retles<double,History,LARGESIZE>(
                n,work_.data(),out,a,acc,b,bcc,integrator_,beta_,dt_);
        else
            cntr::convolution_timestep_les_lesadv<double,History,LARGESIZE>(
                n,work_.data(),out,a,acc,b,bcc,integrator_,beta_,dt_);
    }
    void lesser_product(History& out,History& a,History& acc,
                        History& b,History& bcc,bool retles) {
        for (int n=0;n<=nt_;++n) {
            lesser_slice(n,out,a,acc,b,bcc,retles);
            std::copy_n(work_.data(),size_t(n+1)*out.element_size(),out.lesptr(0,n));
        }
    }
public:
    Solver(const std::vector<double>& times,int dimension,int statistics,int order,double beta)
        : order_(order), dimension_(dimension), statistics_(statistics), nt_(int(times.size()-1)),
          dt_(supplied::uniform_step(times,order)), beta_(beta),
          integrator_(integration::I<double>(order)),
          work_(times.size()*size_t(dimension)*size_t(dimension)) {
        if (dimension<1 || (statistics!=-1 && statistics!=1) || !(beta>0.) || !std::isfinite(beta))
            throw std::invalid_argument("invalid dimension/statistics/inert beta");
    }
    void solve(History& result,History& source,History& kernel) {
        check(result); check(source); check(kernel);
        // F=-(source*kernel)^R and Fcc=-(kernel*source)^R. The retarded
        // VIE consumes neither lesser nor imaginary channels. Destroy both
        // kernels before allocating the two retarded dressing operators.
        {
            History f(nt_,dimension_,statistics_,order_,true,false);
            History fcc(nt_,dimension_,statistics_,order_,true,false);
            ret_product(f,source,source,kernel,kernel);
            ret_product(fcc,kernel,kernel,source,source);
            f.negate_ret(); fcc.negate_ret();
            cntr::vie2_start_ret<double,History,LARGESIZE>(result,f,fcc,source,integrator_,dt_);
            // A54 intentionally uses the reverse kernel as the first argument
            // of the right-propagating timestep solver. Do not interchange it.
            for (int n=order_+1;n<=nt_;++n)
                cntr::vie2_timestep_ret<double,History,LARGESIZE>(n,result,fcc,f,source,integrator_,dt_);
        }
        History left(nt_,dimension_,statistics_,order_,false,true);
        History right(nt_,dimension_,statistics_,order_,false,true);
        {
            History a(nt_,dimension_,statistics_,order_,true,false);
            History acc(nt_,dimension_,statistics_,order_,true,false);
            ret_product(a,result,result,kernel,kernel); // A=R*Sigma^R
            ret_product(acc,kernel,kernel,result,result); // Acc=Sigma^R*R
            lesser_product(left,a,acc,source,source,true);
            lesser_product(right,source,source,acc,a,false);
            // q + A*q + q*Acc + (A*q)*Acc. There is no retarded part of
            // q/left/right. Keep BOTH lesser companion histories until the complete
            // sandwich is accumulated; startup reads future nodes through order.
            for (int n=0;n<=nt_;++n) {
                lesser_slice(n,result,left,right,acc,a,false);
                const size_t count=size_t(n+1)*result.element_size();
                auto* out=result.lesptr(0,n);
                const auto* q=source.lesptr(0,n);
                const auto* l=left.lesptr(0,n);
                const auto* r=right.lesptr(0,n);
                for (size_t z=0;z<count;++z) out[z]=q[z]+l[z]+r[z]+work_[z];
            }
        } // A/Acc are dead before noise transport; release both retarded buffers.
        // Reuse the two lesser buffers for R*Sigma^< and Sigma^<*R^A.
        // Neither real-axis helper reads the lesser component of result here.
        lesser_product(left,result,result,kernel,kernel,true);
        lesser_product(right,kernel,kernel,result,result,false);
        for (int n=0;n<=nt_;++n) {
            lesser_slice(n,result,left,right,result,result,false);
            const size_t count=size_t(n+1)*result.element_size();
            auto* out=result.lesptr(0,n);
            for (size_t z=0;z<count;++z) out[z]+=work_[z];
        }
        // Discarded A54 Diagnostics/residual products are deliberately absent.
        // No clipping, exchange symmetrization, truncation, or state repair.
    }
};
}
#endif
