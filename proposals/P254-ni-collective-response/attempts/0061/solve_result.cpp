// Production transport only. The complete admitted A54 algorithm is unchanged.
#include "../0054/supplied_keldysh.hpp"
#include "../0049/component_io.hpp"
#include <iostream>
#include <limits>
#include <filesystem>
#include <charconv>
int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) throw std::invalid_argument("usage: solve_result INPUT.bin RESULT.bin [MAX_INPUT_BYTES]");
        if (std::filesystem::exists(argv[2])) throw std::invalid_argument("result must be new");
        std::ifstream in(argv[1], std::ios::binary);
        const auto h = p254::supplied::read_header(in);
        // Limits are caller policy, not the pilot's allocation. Check exact
        // extents and arithmetic before any complete-history Green allocation.
        uint64_t cap = std::numeric_limits<uint64_t>::max();
        if (argc == 4) {
            const char* end = argv[3] + std::strlen(argv[3]);
            const auto parsed = std::from_chars(argv[3], end, cap);
            if (parsed.ec != std::errc{} || parsed.ptr != end || cap == 0)
                throw std::invalid_argument("MAX_INPUT_BYTES must be a positive uint64");
        }
        const uint64_t count = uint64_t(h.nt) + 1, d = h.dimension;
        const uint64_t cells = count*(count+1)/2;
        const uint64_t base = 32 + count*8;
        const uint64_t max_extent = std::numeric_limits<std::streamoff>::max();
        if (d > uint64_t(std::numeric_limits<int>::max()) / d ||
            cells > (max_extent-base)/64/d/d)
            throw std::invalid_argument("input layout exceeds native integer representation");
        const uint64_t expected = base + cells*d*d*64;
        if (expected > cap || std::filesystem::file_size(argv[1]) != expected)
            throw std::invalid_argument("input byte extent mismatch or caller byte limit");
        std::vector<double> times(count);
        for (auto& t : times) p254::supplied::read_value(in, t);
        p254::supplied::uniform_step(times, h.order);
        const int ntau = std::max(2*h.order+2, 12);
        p254::Green source(h.nt,ntau,h.dimension,h.statistics), kernel(h.nt,ntau,h.dimension,h.statistics), result(h.nt,ntau,h.dimension,h.statistics);
        p254::clear(source); p254::clear(kernel); p254::clear(result);
        for (auto* g : {&source,&kernel}) for (bool lesser : {false,true})
            p254::supplied::read_component(in,*g,lesser);
        p254::keldysh::Solver solver(times,h.dimension,h.statistics,h.order,h.beta);
        solver.solve(result,source,kernel); // physical omission=0, unchanged math
        std::ofstream out(argv[2],std::ios::binary);
        p254::supplied::write_header(out,h);
        for (double t : times) p254::supplied::write_value(out,t);
        for (bool lesser : {false,true}) p254::supplied::write_component(out,result,lesser);
        out.close(); if (!out) throw std::runtime_error("result close failed");
        std::cout << "A61_RESULT_COMPLETE nt=" << h.nt << " dimension=" << h.dimension
                  << " statistics=" << h.statistics << " bytes=" << std::filesystem::file_size(argv[2]) << "\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 2; }
}
