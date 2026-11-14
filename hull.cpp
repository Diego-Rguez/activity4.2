#include "hull.hpp"

std::vector<Pt> graham(std::vector<Pt> P) {
    // 1. remove duplicates; handle n < 3
    std::sort(P.begin(), P.end());
    P.erase(std::unique(P.begin(), P.end()), P.end());
    if (P.size() < 3) return P;
    // 2. p0 = lowest point
    Pt p0 = P[0];
    // 3. sort the rest by angle around p0 using orient() (no atan2!), ties by distance
    std::sort(P.begin() + 1, P.end(), [p0](const Pt& a, const Pt& b) {
        long long o = orient(p0, a, b);
        if (o == 0) return dist2(p0, a) < dist2(p0, b);
        return o > 0;
    });
    // 4. scan with a stack; pop while the last two points and the new one do NOT make a left turn
    std::vector<Pt> H = {};
    for (const Pt& p : P) {
        while (H.size() >= 2 && orient(H[H.size() - 2], H.back(), p) <= 0) {
            H.pop_back();
        }
        H.push_back(p);
    }
    return H; 
}

std::vector<Pt> jarvis(std::vector<Pt> P, long long* tests) {
    std::sort(P.begin(), P.end());
    P.erase(std::unique(P.begin(), P.end()), P.end());
    if (P.size() < 3) return P;
    // start at the lowest point; from the current point choose the candidate such that
    // every other point is to its LEFT (use orient); on ties keep the farthest point
    Pt start = P[0];
    Pt current = start;
    std::vector<Pt> H = {};
    long long cnt = 0;
    do {
        H.push_back(current);
        Pt next = P[1]; // candidate for the next point in the hull
        for (const Pt& p : P) {
            if (p == current) continue;
            long long o = orient(current, next, p);
            cnt++;
            if (o < 0 || (o == 0 && dist2(current, p) > dist2(current, next))) {
                next = p;
            }
        }
        current = next;
    } while (!(current == start));   
    if (tests) *tests = cnt;        
    return H;
}
