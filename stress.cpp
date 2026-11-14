// Part 3: compare your implementations with CGAL on thousands of random inputs.
#include "hull.hpp"
#include "gen.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <iostream>
using K  = CGAL::Exact_predicates_inexact_constructions_kernel;
using P2 = K::Point_2;

std::vector<Pt> cgal_hull(const std::vector<Pt>& P) {
    std::vector<P2> in, out;
    for (auto& p : P) in.emplace_back((double)p.x, (double)p.y);   // exact: |coords| <= 1e9 < 2^53
    CGAL::convex_hull_2(in.begin(), in.end(), std::back_inserter(out));
    std::vector<Pt> H;
    for (auto& q : out) H.push_back({(long long)q.x(), (long long)q.y()});
    // rotate H so that it starts at the lowest point (our output convention)
    std::rotate(H.begin(), std::min_element(H.begin(), H.end()), H.end());
    return H;
}

int main() {
    int failures = 0;
    for (unsigned seed = 0; seed < 1000; seed++)
        for (std::string d : {"square", "disk", "parabola"}) {
            int n = 1 + seed % 300;
            auto P = gen(d, n, seed);
            // Hard cases for the stress test: duplicated points and collinear points on an edge
            if (seed % 3 == 1) {
                size_t k = std::min<size_t>(5, P.size());
                for (size_t i = 0; i < k; i++) P.push_back(P[i]);
            } else if (seed % 3 == 2) {
                for (long long k = 0; k <= 4; k++) P.push_back({-1000 + 500 * k, -1000000000LL});
            }
            auto C = cgal_hull(P);
            if (graham(P) != C || jarvis(P) != C) {
                if (++failures <= 5) std::cout << "mismatch: dist=" << d << " n=" << n << " seed=" << seed << "\n";
            }
        }
    std::cout << (failures ? "FAILED: " : "OK: ") << failures << " mismatches\n";
}
