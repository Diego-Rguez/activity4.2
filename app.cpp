// Part 5: "Fence the orchard". Input format as always. Output (3 decimals):
//   h, perimeter, area, diameter (largest distance between two trees), and the pair attaining it.
// "./app brute" checks the rotating calipers against brute force and times both (question 5b).
#include "hull.hpp"
#include "gen.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// Area with the shoelace formula. We return TWICE the area so everything stays
// in integers and the result is exact (the /2 only happens when printing).
long long twice_area(const std::vector<Pt>& H) {
    long long sum = 0;
    int h = H.size();
    for (int i = 0; i < h; i++) {
        const Pt& p = H[i];
        const Pt& q = H[(i + 1) % h];   // next vertex; the last one wraps around to the first
        sum += p.x * q.y - q.x * p.y;
    }
    return sum;   // positive because the hull comes in counter-clockwise order
}

// Perimeter: the length of every hull edge added up, including the one that closes the polygon.
double perimeter(const std::vector<Pt>& H) {
    double total = 0;
    int h = H.size();
    for (int i = 0; i < h; i++) {
        total += std::sqrt((double)dist2(H[i], H[(i + 1) % h]));
    }
    return total;
}

// Squared diameter with rotating calipers in O(h). The farthest pair is stored in a, b.
long long diameter2(const std::vector<Pt>& H, Pt& a, Pt& b) {
    int h = H.size();
    if (h == 0) return 0;
    long long best = -1;
    int j = 1 % h;   // second pointer (it's 0 when there is only one point)
    for (int i = 0; i < h; i++) {
        int ni = (i + 1) % h;
        // move j forward while the next vertex is farther from the edge H[i] -> H[ni]
        while (orient(H[i], H[ni], H[(j + 1) % h]) > orient(H[i], H[ni], H[j])) {
            j = (j + 1) % h;
        }
        // H[j] is the farthest vertex from this edge, so check it against both endpoints
        if (dist2(H[i], H[j]) > best)  { best = dist2(H[i], H[j]);  a = H[i];  b = H[j]; }
        if (dist2(H[ni], H[j]) > best) { best = dist2(H[ni], H[j]); a = H[ni]; b = H[j]; }
    }
    return best;
}

// Brute force diameter: try every pair of points, O(n^2). Only used to check diameter2.
long long brute_diameter2(const std::vector<Pt>& P) {
    long long best = 0;
    for (size_t i = 0; i < P.size(); i++)
        for (size_t j = i + 1; j < P.size(); j++)
            best = std::max(best, dist2(P[i], P[j]));
    return best;
}

double ms_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

int brute_mode() {
    // 1) Validation: on random inputs both methods must give the same squared diameter.
    //    We compare the distance and not the pair, because two different pairs can tie.
    const std::string dists[] = {"square", "disk", "parabola"};
    int mismatches = 0;
    for (int t = 0; t < 3000; t++) {
        auto P = gen(dists[t % 3], 1 + t % 200, t);
        Pt a{0, 0}, b{0, 0};
        if (diameter2(graham(P), a, b) != brute_diameter2(P)) mismatches++;
    }
    std::cout << "validation: " << mismatches << " mismatches in 3000 random inputs\n";

    // 2) Timing with n = 1e5 points in a square.
    auto P = gen("square", 100000, 2038);
    auto start = std::chrono::steady_clock::now();
    long long d_brute = brute_diameter2(P);
    double t_brute = ms_since(start);

    start = std::chrono::steady_clock::now();
    Pt a{0, 0}, b{0, 0};
    auto H = graham(P);
    long long d_calipers = diameter2(H, a, b);
    double t_calipers = ms_since(start);

    std::cout << std::fixed << std::setprecision(3)
              << "n = " << P.size() << ", h = " << H.size() << "\n"
              << "brute force O(n^2):        " << t_brute << " ms\n"
              << "hull + calipers O(n log n): " << t_calipers << " ms\n"
              << (d_brute == d_calipers ? "same diameter" : "DIFFERENT diameter!") << "\n";
    return mismatches == 0 && d_brute == d_calipers ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "brute") return brute_mode();

    int n; std::cin >> n;
    std::vector<Pt> P(n);
    for (auto& p : P) std::cin >> p.x >> p.y;
    auto H = graham(P);
    Pt a{0, 0}, b{0, 0};
    long long d2 = diameter2(H, a, b);
    std::cout << std::fixed << std::setprecision(3)
              << "h = " << H.size() << "\nperimeter = " << perimeter(H)
              << "\narea = " << twice_area(H) / 2.0 << "\ndiameter = " << std::sqrt((double)d2)
              << "  between (" << a.x << ", " << a.y << ") and (" << b.x << ", " << b.y << ")\n";
}
