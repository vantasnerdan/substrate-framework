// P254VIE1 result-only executable; real-axis ABI is identical to A61.
#include "real_axis_keldysh.hpp"
#include <iostream>
#include <filesystem>
#include <charconv>
int main(int argc,char** argv) {
    try {
        if (argc!=3 && argc!=4)
            throw std::invalid_argument("usage: solve_result INPUT.bin RESULT.bin [MAX_INPUT_BYTES]");
        if (std::filesystem::exists(argv[2])) throw std::invalid_argument("result must be new");
        std::ifstream in(argv[1],std::ios::binary);
        const auto h=p254::supplied::read_header(in);
        uint64_t cap=std::numeric_limits<uint64_t>::max();
        if (argc==4) {
            const char* end=argv[3]+std::strlen(argv[3]);
            const auto parsed=std::from_chars(argv[3],end,cap);
            if (parsed.ec!=std::errc{} || parsed.ptr!=end || cap==0)
                throw std::invalid_argument("MAX_INPUT_BYTES must be a positive uint64");
        }
        const uint64_t count=uint64_t(h.nt)+1,d=h.dimension;
        const uint64_t cells=count*(count+1)/2,base=32+count*8;
        const uint64_t max_extent=std::numeric_limits<std::streamoff>::max();
        if (d>uint64_t(std::numeric_limits<int>::max())/d ||
            cells>(max_extent-base)/64/d/d)
            throw std::invalid_argument("input layout exceeds native integer representation");
        const uint64_t expected=base+cells*d*d*64;
        if (expected>cap || std::filesystem::file_size(argv[1])!=expected)
            throw std::invalid_argument("input byte extent mismatch or caller byte limit");
        std::vector<double> times(count);
        for (auto& t:times) p254::supplied::read_value(in,t);
        p254::supplied::uniform_step(times,h.order);
        // Only R and < can be imported by this schema. No imaginary branch is
        // allocated, inferred, or silently ignored by the native representation.
        p254::real_axis::History source(h.nt,h.dimension,h.statistics,h.order,true,true);
        p254::real_axis::History kernel(h.nt,h.dimension,h.statistics,h.order,true,true);
        p254::real_axis::History result(h.nt,h.dimension,h.statistics,h.order,true,true);
        for (auto* g:{&source,&kernel}) for (bool lesser:{false,true})
            p254::real_axis::read_component(in,*g,lesser);
        p254::real_axis::Solver solver(times,h.dimension,h.statistics,h.order,h.beta);
        solver.solve(result,source,kernel);
        std::ofstream out(argv[2],std::ios::binary);
        p254::supplied::write_header(out,h);
        for (double t:times) p254::supplied::write_value(out,t);
        for (bool lesser:{false,true}) p254::real_axis::write_component(out,result,lesser);
        out.close(); if (!out) throw std::runtime_error("result close failed");
        std::cout<<"A85_RESULT_COMPLETE nt="<<h.nt<<" dimension="<<h.dimension
                 <<" statistics="<<h.statistics<<" bytes="<<std::filesystem::file_size(argv[2])<<"\n";
        return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<"\n"; return 2; }
}
