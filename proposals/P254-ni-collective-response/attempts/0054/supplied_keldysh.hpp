// Fixed supplied-kernel Dyson: upstream retarded VIE plus two-sided Langreth lesser.
#ifndef P254_SUPPLIED_KELDYSH_HPP
#define P254_SUPPLIED_KELDYSH_HPP
#include "../0049/supplied_vie2.hpp"

namespace p254::keldysh {
struct Diagnostics {
    double lesser_exchange=0., retarded_residual=0., lesser_left_residual=0., lesser_right_residual=0.;
    double initial_term=0., noise_term=0., product_order_difference=0.;
};
class Solver {
    int order_;
    double dt_, beta_;
    supplied::Solver kernels_;
    Green make(const Green& g) const {
        Green out(g.nt(),g.ntau(),g.size1(),g.sig()); clear(out); return out;
    }
    void product(Green& c, Green& a, Green& acc, Green& b, Green& bcc) const {
        for (int n=0;n<=c.nt();++n)
            cntr::convolution_timestep(n,c,a,acc,b,bcc,integration::I<double>(order_),beta_,dt_);
    }
    static void component(Green& dest, const Green& source, bool lesser) {
        const int size=source.element_size();
        for (int n=0;n<=source.nt();++n) for (int m=0;m<=n;++m) {
            const auto* from=lesser ? source.lesptr(m,n) : source.retptr(n,m);
            auto* to=lesser ? dest.lesptr(m,n) : dest.retptr(n,m);
            std::copy(from,from+size,to);
        }
    }
public:
    Solver(const std::vector<double>& times,int dimension,int statistics,int order,double beta=1.)
        : order_(order),dt_(supplied::uniform_step(times,order)),beta_(beta),
          kernels_(times,dimension,statistics,order,beta) {}
    // Full history, including 0..order, must be imported. No selfconsistent update.
    // omit: 0 physical; 1 omit initial-source; 2 omit noise; 3 swap A/Acc orders.
    Diagnostics solve(Green& result,Green& source,Green& kernel,int omit=0) {
        kernels_.check(result); kernels_.check(source); kernels_.check(kernel);
        if (omit<0 || omit>3) throw std::invalid_argument("unknown explicit diagnostic omission");
        for (int n=0;n<=result.nt();++n) kernels_.refresh(n,source,kernel);
        auto& f=const_cast<Green&>(kernels_.forward_kernel());
        auto& fcc=const_cast<Green&>(kernels_.reverse_kernel());
        auto& i=integration::I<double>(order_);
        cntr::vie2_start_ret<double,Green,LARGESIZE>(result,f,fcc,source,i,dt_);
        for (int n=order_+1;n<=result.nt();++n)
            cntr::vie2_timestep_ret<double,Green,LARGESIZE>(n,result,fcc,f,source,i,dt_);
        Green r=make(source),s=make(source),q=make(source),k=make(source);
        component(r,result,false); component(s,kernel,false);
        component(q,source,true); component(k,kernel,true);
        Green a=make(source),acc=make(source),left=make(source),right=make(source);
        Green sandwich=make(source),nk=make(source),nkcc=make(source),noise=make(source);
        product(a,r,r,s,s); product(acc,s,s,r,r); // genuinely independent orders
        auto& al=(omit==3 ? acc : a);
        auto& ac=(omit==3 ? a : acc);
        product(left,al,ac,q,q); product(right,q,q,ac,al);
        product(sandwich,left,right,ac,al);
        product(nk,r,r,k,k); product(nkcc,k,k,r,r);
        product(noise,nk,nkcc,r,r);
        Diagnostics d;
        Matrix x(source.size1(),source.size1());
        for (int n=0;n<=result.nt();++n) for (int m=0;m<=n;++m) {
            auto* out=result.lesptr(m,n);
            for (int z=0;z<source.element_size();++z) {
                const auto initial=q.lesptr(m,n)[z]+left.lesptr(m,n)[z]+right.lesptr(m,n)[z]+sandwich.lesptr(m,n)[z];
                d.initial_term=std::max(d.initial_term,std::abs(initial));
                d.noise_term=std::max(d.noise_term,std::abs(noise.lesptr(m,n)[z]));
                out[z]=(omit==1 ? Complex(0.) : initial)+(omit==2 ? Complex(0.) : noise.lesptr(m,n)[z]);
                d.product_order_difference=std::max(d.product_order_difference,std::abs(a.retptr(n,m)[z]-acc.retptr(n,m)[z]));
            }
            if (n==m) {
                result.get_les(n,n,x);
                d.lesser_exchange=std::max(d.lesser_exchange,(x+x.adjoint()).cwiseAbs().maxCoeff());
            }
        }
        // Raw finite-step VIE residuals are measurements, never a repair target.
        Green residual=make(source),residualcc=make(source);
        product(residual,f,fcc,result,result);
        product(residualcc,result,result,fcc,f);
        for (int n=0;n<=result.nt();++n) for (int m=0;m<=n;++m)
            for (int z=0;z<source.element_size();++z) {
                d.retarded_residual=std::max(d.retarded_residual,std::abs(result.retptr(n,m)[z]+residual.retptr(n,m)[z]-source.retptr(n,m)[z]));
                d.lesser_left_residual=std::max(d.lesser_left_residual,std::abs(result.lesptr(m,n)[z]+residual.lesptr(m,n)[z]-source.lesptr(m,n)[z]));
                d.lesser_right_residual=std::max(d.lesser_right_residual,std::abs(result.lesptr(m,n)[z]+residualcc.lesptr(m,n)[z]-source.lesptr(m,n)[z]));
            }
        return d;
    }
    const Green& forward_kernel() const { return kernels_.forward_kernel(); }
    const Green& reverse_kernel() const { return kernels_.reverse_kernel(); }
};
}
#endif
