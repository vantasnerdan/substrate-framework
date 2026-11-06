// Persist supplied Dyson result and raw finite-step diagnostics, never repair data.
#include "supplied_keldysh.hpp"
#include "../0049/component_io.hpp"
#include <iostream>
#include <iomanip>
int main(int argc,char** argv) {
    try {
        if (argc!=3 && argc!=4) throw std::invalid_argument("usage: solve_components INPUT.bin OUTPUT.bin [EXPLICIT_OMISSION_1_TO_3]");
        const int omit=argc==4 ? std::stoi(argv[3]) : 0;
        std::ifstream in(argv[1],std::ios::binary);
        auto h=p254::supplied::read_header(in);
        std::vector<double> times(h.nt+1);
        for (auto& t:times) p254::supplied::read_value(in,t);
        const int ntau=std::max(2*h.order+2,12);
        p254::Green source(h.nt,ntau,h.dimension,h.statistics),kernel(h.nt,ntau,h.dimension,h.statistics),result(h.nt,ntau,h.dimension,h.statistics);
        p254::clear(source); p254::clear(kernel); p254::clear(result);
        for (auto* g:{&source,&kernel}) for (bool lesser:{false,true}) p254::supplied::read_component(in,*g,lesser);
        if (in.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("trailing input bytes");
        p254::keldysh::Solver solver(times,h.dimension,h.statistics,h.order,h.beta);
        auto d=solver.solve(result,source,kernel,omit);
        std::ofstream out(argv[2],std::ios::binary);
        p254::supplied::write_header(out,h);
        for (double t:times) p254::supplied::write_value(out,t);
        for (const auto* g:std::initializer_list<const p254::Green*>{&result,&solver.forward_kernel(),&solver.reverse_kernel()})
            for (bool lesser:{false,true}) p254::supplied::write_component(out,*g,lesser);
        out.close(); if (!out) throw std::runtime_error("output close failed");
        std::ofstream diag(std::string(argv[2])+".diagnostics.json");
        diag<<std::setprecision(17)<<"{\n  \"omission\": "<<omit
            <<",\n  \"lesser_exchange_max\": "<<d.lesser_exchange
            <<",\n  \"retarded_vie_residual_max\": "<<d.retarded_residual
            <<",\n  \"lesser_left_vie_residual_max\": "<<d.lesser_left_residual
            <<",\n  \"lesser_right_vie_residual_max\": "<<d.lesser_right_residual
            <<",\n  \"initial_source_term_max\": "<<d.initial_term
            <<",\n  \"noise_term_max\": "<<d.noise_term
            <<",\n  \"retarded_product_order_difference_max\": "<<d.product_order_difference<<"\n}\n";
        diag.close(); if (!diag) throw std::runtime_error("diagnostic close failed");
        std::cout<<"supplied_keldysh_complete nt="<<h.nt<<" dimension="<<h.dimension<<" statistics="<<h.statistics<<" omission="<<omit<<"\n";
        return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<"\n"; return 2; }
}
