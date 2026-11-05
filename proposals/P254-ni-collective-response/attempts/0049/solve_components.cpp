// Actual supplied-kernel consumer, intentionally separate from Gaussian producers.
#include "component_io.hpp"
#include <iostream>
int main(int argc,char** argv) {
    try {
        if (argc!=3) throw std::invalid_argument("usage: solve_components INPUT.bin OUTPUT.bin");
        std::ifstream in(argv[1],std::ios::binary);
        auto h=p254::supplied::read_header(in);
        std::vector<double> times(h.nt+1);
        for (auto& t:times) p254::supplied::read_value(in,t);
        p254::supplied::uniform_step(times,h.order);
        const int ntau=std::max(2*h.order+2,12);
        p254::Green source(h.nt,ntau,h.dimension,h.statistics),kernel(h.nt,ntau,h.dimension,h.statistics),result(h.nt,ntau,h.dimension,h.statistics);
        p254::clear(source); p254::clear(kernel); p254::clear(result);
        for (auto* g:{&source,&kernel}) for (bool lesser:{false,true})
            p254::supplied::read_component(in,*g,lesser);
        if (in.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("trailing input bytes");
        p254::supplied::Solver solver(times,h.dimension,h.statistics,h.order,h.beta);
        solver.start(result,source,kernel);
        for (int n=h.order+1;n<=h.nt;++n) solver.step(n,result,source,kernel);
        std::ofstream out(argv[2],std::ios::binary);
        p254::supplied::write_header(out,h);
        for (double t:times) p254::supplied::write_value(out,t);
        for (const auto* g:std::initializer_list<const p254::Green*>{&result,&solver.forward_kernel(),&solver.reverse_kernel()})
            for (bool lesser:{false,true}) p254::supplied::write_component(out,*g,lesser);
        out.close(); if (!out) throw std::runtime_error("output close failed");
        std::cout<<"supplied_vie2_complete nt="<<h.nt<<" dimension="<<h.dimension<<" statistics="<<h.statistics<<"\n";
        return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<"\n"; return 2; }
}
