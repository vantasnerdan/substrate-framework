// Real-time, two-branch Gaussian IVP using pinned NESSi's unmodified VIE2.
// Not an L-contour equilibrium preparation: g^< carries the original rho.
#ifndef P254_NESSI_GAUSSIAN_IVP_HPP
#define P254_NESSI_GAUSSIAN_IVP_HPP
#include <cntr/cntr.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <complex>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace p254 {
using Matrix = Eigen::MatrixXcd;
using Green = cntr::herm_matrix<double>;
using Complex = std::complex<double>;
const Complex imaginary(0.0, 1.0);

struct GaussianInput {
    int nt, order;
    double dt;
    Matrix rho_system, rho_bath;
    // H[n] is held on [n*dt,(n+1)*dt); V[n] is its nodal vertex.
    // Last H node is retained for provenance, with no subsequent interval.
    std::vector<Matrix> h_system, h_bath, coupling;
};

inline void clear(Green &g) {
    for (int n = -1; n <= g.nt(); ++n) g.set_timestep_zero(n);
}

inline void require_hermitian(const Matrix &a) {
    if (!a.allFinite() || (a-a.adjoint()).norm() > 1.e-10*(1.0+a.norm()))
        throw std::invalid_argument("nonfinite/non-Hermitian Hamiltonian or state; no repair");
}

inline void validate(const GaussianInput &x) {
    const int ns = x.rho_system.rows(), nb = x.rho_bath.rows();
    if (ns < 1 || nb < 1 || x.rho_system.cols()!=ns || x.rho_bath.cols()!=nb ||
        x.order < 1 || x.order > 5 || x.nt < x.order || !(x.dt>0) || !std::isfinite(x.dt))
        throw std::invalid_argument("invalid IVP dimensions, uniform grid or order (1..5)");
    if (x.h_system.size()!=size_t(x.nt+1) || x.h_bath.size()!=size_t(x.nt+1) ||
        x.coupling.size()!=size_t(x.nt+1)) throw std::invalid_argument("missing time nodes");
    for (const Matrix *rho : {&x.rho_system, &x.rho_bath}) {
        require_hermitian(*rho);
        Eigen::SelfAdjointEigenSolver<Matrix> eig(*rho, Eigen::EigenvaluesOnly);
        if (eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff() < -1.e-10 ||
            eig.eigenvalues().maxCoeff() > 1.+1.e-10)
            throw std::invalid_argument("not a fermionic Gaussian density; no clipping");
    }
    for (int n=0;n<=x.nt;++n) {
        if (x.h_system[n].rows()!=ns || x.h_system[n].cols()!=ns ||
            x.h_bath[n].rows()!=nb || x.h_bath[n].cols()!=nb ||
            x.coupling[n].rows()!=ns || x.coupling[n].cols()!=nb ||
            !x.coupling[n].allFinite()) throw std::invalid_argument("invalid operator shape");
        require_hermitian(x.h_system[n]);
        require_hermitian(x.h_bath[n]);
    }
}

inline std::vector<Matrix> free_evolution(const std::vector<Matrix> &h, double dt) {
    std::vector<Matrix> u(h.size());
    u[0]=Matrix::Identity(h[0].rows(),h[0].rows());
    for (size_t n=1;n<h.size();++n) u[n]=(-imaginary*dt*h[n-1]).exp()*u[n-1];
    return u;
}

inline void fill_free(Green &g, const std::vector<Matrix> &u, const Matrix &rho) {
    const int ns=rho.rows();
    Matrix ret(ns,ns), les(ns,ns), occupied_right(ns,ns);
    for (int n=0;n<=g.nt();++n) {
        occupied_right.noalias()=rho*u[n].adjoint();
        for (int m=0;m<=n;++m) {
            ret.noalias()=-imaginary*u[n]*u[m].adjoint();
            les.noalias()=imaginary*u[m]*occupied_right;
            g.set_ret(n,m,ret);
            g.set_les(m,n,les);
        }
    }
}

inline void fill_embedding(Green &sigma, const std::vector<Matrix> &u,
                           const Matrix &rho, const std::vector<Matrix> &v) {
    // Sigma(t,s)=V(t) g_b(t,s) V(s)^dagger, including lesser injection.
    std::vector<Matrix> vu(u.size());
    for (size_t n=0;n<u.size();++n) vu[n].noalias()=v[n]*u[n];
    const int ns=v.front().rows(), nb=rho.rows();
    Matrix ret(ns,ns), les(ns,ns), occupied_right(nb,ns);
    for (int n=0;n<=sigma.nt();++n) {
        occupied_right.noalias()=rho*vu[n].adjoint();
        for (int m=0;m<=n;++m) {
            ret.noalias()=-imaginary*vu[n]*vu[m].adjoint();
            les.noalias()=imaginary*vu[m]*occupied_right;
            sigma.set_ret(n,m,ret);
            sigma.set_les(m,n,les);
        }
    }
}

inline void solve_real_time_vie2(Green &g, Green &g0, Green &sigma,
                                int order, double dt, double inert_beta=1.0) {
    // On the real-time closed contour, G=g0+g0*Sigma*G. Its Langreth
    // products contain only ret*ret and ret*less+less*adv integrals from 0.
    // Initial product Gaussian correlations are ALREADY in g0^< and g_b^<.
    // The unused imaginary storage is NOT an equilibrium contour function.
    // All tv legs are zero; upstream convolution tv=Atv*Bmat+Aret*Btv,
    // and the only extra lesser term is Atv*Bvt (convolution_impl 901,1076).
    // Hence zero tv is an invariant subspace, regardless of mat/beta.
    // DO NOT use dyson_start: dyson_timestep_les initializes from Gtv(n,0).
    const int nt=g.nt(), ntau=g.ntau(), ns=g.size1();
    Green f(nt,ntau,ns,-1), fcc(nt,ntau,ns,-1);
    clear(f); clear(fcc);
    auto &integrator=integration::I<double>(order);
    for (int n=0;n<=nt;++n) {
        cntr::convolution_timestep(n,f,g0,sigma,integrator,inert_beta,dt);
        cntr::convolution_timestep(n,fcc,sigma,g0,integrator,inert_beta,dt);
        f.smul(n,-1.0); fcc.smul(n,-1.0);
    }
    // Explicit startup bypasses vie2(), which would solve a Matsubara problem.
    cntr::vie2_start(g,f,fcc,g0,integrator,inert_beta,dt);
    for (int n=order+1;n<=nt;++n)
        cntr::vie2_timestep(n,g,f,fcc,g0,integrator,inert_beta,dt);
}

inline void solve_embedding(const GaussianInput &x, Green &g, Green &sigma,
                            double inert_beta=1., double inert_mat=0.) {
    validate(x);
    if (!(inert_beta>0.) || !std::isfinite(inert_beta) || !std::isfinite(inert_mat))
        throw std::invalid_argument("invalid inert-storage probe");
    const int ntau=std::max(2*x.order+2,12), ns=x.rho_system.rows();
    g.resize_discard(x.nt,ntau,ns); sigma.resize_discard(x.nt,ntau,ns);
    g.set_sig(-1); sigma.set_sig(-1);
    Green g0(x.nt,ntau,ns,-1);
    clear(g); clear(g0); clear(sigma);
    fill_free(g0,free_evolution(x.h_system,x.dt),x.rho_system);
    fill_embedding(sigma,free_evolution(x.h_bath,x.dt),x.rho_bath,x.coupling);
    // Exposing probe: arbitrary imaginary arrays, never a thermal preparation.
    Matrix sentinel=inert_mat*Matrix::Identity(ns,ns);
    for (int tau=0;tau<=ntau;++tau) {
        g.set_mat(tau,sentinel); g0.set_mat(tau,sentinel); sigma.set_mat(tau,sentinel);
    }
    solve_real_time_vie2(g,g0,sigma,x.order,x.dt,inert_beta);
}
}
#endif
