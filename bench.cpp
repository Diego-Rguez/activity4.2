// Part 4: timing experiments. Prints CSV: dist,n,h,algo,ms,tests
#include "hull.hpp"
#include "gen.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_jarvis.h>
#include <chrono>
#include <iostream>
using K  = CGAL::Exact_predicates_inexact_constructions_kernel;
using P2 = K::Point_2;

template <class F> double time_ms(F f) {
    auto t0 = std::chrono::steady_clock::now(); f();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

int main() {
    std::cout << "dist,n,h,algo,ms,tests\n";
    for (std::string d : {"square", "disk", "parabola"})
        for (int n : {1000, 3000, 10000, 30000, 100000, 300000, 1000000}) {
            if (d == "parabola" && n > 60000) continue;
            auto P = gen(d, n, 2026);
            std::vector<P2> in;
            for (auto& p : P) in.emplace_back((double)p.x, (double)p.y);
            size_t h = 0;
            double t = time_ms([&] { h = graham(P).size(); });
            std::cout << d << "," << n << "," << h << ",graham," << t << ",\n";
            // Jarvis costs about n*h steps, so we skip both versions when that gets too big
            bool run_jarvis = (long long)n * h <= 2000000000LL;
            if (run_jarvis) {
                long long tests = 0;
                t = time_ms([&] { jarvis(P, &tests); });
                std::cout << d << "," << n << "," << h << ",jarvis," << t << "," << tests << "\n";
            }

            // CGAL versions, each one writes its hull into a fresh vector
            t = time_ms([&] { std::vector<P2> out; CGAL::convex_hull_2(in.begin(), in.end(), std::back_inserter(out)); });
            std::cout << d << "," << n << "," << h << ",cgal_default," << t << ",\n";

            t = time_ms([&] { std::vector<P2> out; CGAL::ch_graham_andrew(in.begin(), in.end(), std::back_inserter(out)); });
            std::cout << d << "," << n << "," << h << ",cgal_graham_andrew," << t << ",\n";

            if (run_jarvis) {
                t = time_ms([&] { std::vector<P2> out; CGAL::ch_jarvis(in.begin(), in.end(), std::back_inserter(out)); });
                std::cout << d << "," << n << "," << h << ",cgal_jarvis," << t << ",\n";
            }
        }
}
