// Versioned little-endian binary component transport (complex128 row-major).
#ifndef P254_COMPONENT_IO_HPP
#define P254_COMPONENT_IO_HPP
#include "supplied_vie2.hpp"
#include <fstream>
#include <cstdint>
#include <cstring>

namespace p254::supplied {
struct Header { int32_t nt, dimension, order, statistics; double beta; };
template<class T> void read_value(std::istream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value),sizeof(T));
    if (!in) throw std::runtime_error("truncated binary components");
}
template<class T> void write_value(std::ostream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value),sizeof(T));
    if (!out) throw std::runtime_error("binary output failure");
}
inline Header read_header(std::istream& in) {
    uint16_t endian=1;
    if (*reinterpret_cast<char*>(&endian)!=1) throw std::runtime_error("little endian host required");
    char magic[8]; in.read(magic,8);
    if (!in || std::memcmp(magic,"P254VIE1",8)) throw std::runtime_error("unsupported input schema");
    Header h; read_value(in,h.nt); read_value(in,h.dimension);
    read_value(in,h.order); read_value(in,h.statistics); read_value(in,h.beta);
    if (h.nt<h.order || h.order<1 || h.order>5 || h.dimension<1 ||
        (h.statistics!=-1 && h.statistics!=1) || !(h.beta>0.) || !std::isfinite(h.beta))
        throw std::runtime_error("invalid component header");
    return h;
}
inline void write_header(std::ostream& out, const Header& h) {
    out.write("P254VIE1",8);
    write_value(out,h.nt); write_value(out,h.dimension); write_value(out,h.order);
    write_value(out,h.statistics); write_value(out,h.beta);
}
inline void read_component(std::istream& in, Green& g, bool lesser) {
    Matrix a(g.size1(),g.size1());
    for (int n=0;n<=g.nt();++n) for (int m=0;m<=n;++m) {
        for (int r=0;r<g.size1();++r) for (int s=0;s<g.size1();++s) {
            double re,im; read_value(in,re); read_value(in,im); a(r,s)={re,im};
        }
        if (!a.allFinite()) throw std::runtime_error("nonfinite supplied component");
        if (lesser && n==m && (a+a.adjoint()).norm()>1.e-10*(1.+a.norm()))
            throw std::runtime_error("lesser diagonal not contour-Hermitian; no repair");
        if (lesser) g.set_les(m,n,a); else g.set_ret(n,m,a);
    }
}
inline void write_component(std::ostream& out, const Green& g, bool lesser) {
    Matrix a(g.size1(),g.size1());
    for (int n=0;n<=g.nt();++n) for (int m=0;m<=n;++m) {
        if (lesser) g.get_les(m,n,a); else g.get_ret(n,m,a);
        if (!a.allFinite()) throw std::runtime_error("nonfinite solver output");
        for (int r=0;r<g.size1();++r) for (int s=0;s<g.size1();++s) {
            write_value(out,a(r,s).real()); write_value(out,a(r,s).imag());
        }
    }
}
}
#endif
