// Native executable consumer: reads exact Gaussian/operator nodes, emits NESSi
// retarded/lesser/embedding/density arrays. No full-system reference solver here.
#include "nessi_gaussian_ivp.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {
p254::Matrix read_matrix(std::istream &stream, int rows, int cols) {
    p254::Matrix result(rows,cols);
    for (int i=0;i<rows;++i) for (int j=0;j<cols;++j) {
        double real,imag;
        if (!(stream>>real>>imag)) throw std::runtime_error("truncated matrix input");
        result(i,j)={real,imag};
    }
    return result;
}

p254::GaussianInput read_input(const char *path) {
    std::ifstream stream(path);
    std::string magic;
    p254::GaussianInput x;
    int ns,nb;
    if (!(stream>>magic>>ns>>nb>>x.nt>>x.order>>x.dt) ||
        magic!="P254_GAUSSIAN_IVP_V1" || ns<1 || nb<1 || x.nt<1)
        throw std::runtime_error("invalid Gaussian IVP header");
    x.h_system.reserve(x.nt+1);
    x.h_bath.reserve(x.nt+1);
    x.coupling.reserve(x.nt+1);
    x.rho_system=read_matrix(stream,ns,ns);
    x.rho_bath=read_matrix(stream,nb,nb);
    for (int n=0;n<=x.nt;++n) {
        x.h_system.push_back(read_matrix(stream,ns,ns));
        x.h_bath.push_back(read_matrix(stream,nb,nb));
        x.coupling.push_back(read_matrix(stream,ns,nb));
    }
    std::string extra;
    if (stream>>extra) throw std::runtime_error("trailing IVP input");
    return x;
}

void emit_matrix(std::ostream &stream, const p254::Matrix &matrix) {
    for (int i=0;i<matrix.rows();++i) for (int j=0;j<matrix.cols();++j)
        stream<<' '<<matrix(i,j).real()<<' '<<matrix(i,j).imag();
    stream<<'\n';
}

void emit_triangle(std::ostream &stream, const char *name, p254::Green &g, bool lesser) {
    stream<<name<<'\n';
    p254::Matrix matrix(g.size1(),g.size1());
    for (int n=0;n<=g.nt();++n) for (int m=0;m<=n;++m) {
        if (lesser) g.get_les(m,n,matrix); else g.get_ret(n,m,matrix);
        stream<<n<<' '<<m;
        emit_matrix(stream,matrix);
    }
}

void emit_output(const char *path, const p254::GaussianInput &x,
                 p254::Green &g, p254::Green &sigma, double beta, double mat) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot create native output");
    stream<<std::setprecision(std::numeric_limits<double>::max_digits10);
    stream<<"P254_NESSI_REALTIME_V1 "<<g.size1()<<' '<<x.nt<<' '<<x.dt<<'\n';
    emit_triangle(stream,"retarded_n_m",g,false);
    emit_triangle(stream,"lesser_m_n",g,true);
    emit_triangle(stream,"sigma_retarded_n_m",sigma,false);
    emit_triangle(stream,"sigma_lesser_m_n",sigma,true);
    stream<<"density_n"<<'\n';
    p254::Matrix matrix(g.size1(),g.size1());
    for (int n=0;n<=x.nt;++n) {
        g.density_matrix(n,matrix);
        stream<<n;
        emit_matrix(stream,matrix);
    }
    double mixed_max=0.;
    for (int n=0;n<=x.nt;++n) for (int tau=0;tau<=g.ntau();++tau) {
        g.get_tv(n,tau,matrix);
        mixed_max=std::max(mixed_max,matrix.cwiseAbs().maxCoeff());
    }
    stream<<"inactive_mixed_max "<<mixed_max<<'\n';
    stream<<"inert_branch_probe "<<beta<<' '<<mat<<'\n';
    stream.close();
    if (!stream) throw std::runtime_error("failed writing native output");
}
}

int main(int argc,char **argv) {
    try {
        if (argc!=3 && argc!=5)
            throw std::invalid_argument("usage: nessi_gaussian_embedding INPUT OUTPUT [INERT_BETA INERT_MAT]");
        const auto x=read_input(argv[1]);
        const double beta=argc==5?std::stod(argv[3]):1.;
        const double mat=argc==5?std::stod(argv[4]):0.;
        p254::Green g,sigma;
        p254::solve_embedding(x,g,sigma,beta,mat);
        emit_output(argv[2],x,g,sigma,beta,mat);
        std::cout<<"NESSi VIE2 complete: "<<x.nt+1<<" full-memory time nodes; "
                 <<g.size1()<<" system orbitals; original Gaussian rho retained\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr<<"Gaussian IVP error: "<<error.what()<<'\n';
        return 1;
    }
}
