# Lab 4.3 · Convex Hulls in Practice — Report

**TC2038 Analysis and Design of Advanced Algorithms**
Team: Diego Rodríguez Sánchez (A01742085) · Omar Avila Meza (A01743547)

> This file covers Parts 0, 1, 2 and 3. Parts 4, 5 and 6 are added in their own sections.

---

## Part 0 · Setup

CGAL was installed through MSYS2 (`pacman -S mingw-w64-ucrt-x86_64-cgal`). Output of `cgal_hello` on `example13.txt`:

```
5
192 146
312 56
684 86
720 170
420 308
```

These are the 5 expected vertices in counter-clockwise order, but CGAL starts at `(192,146)` instead of the lowest point `(312,56)`.

### Question 0
To compare two implementations automatically, their outputs must be **identical element by element**. The same convex polygon can be written in many ways: it can start at any of its h vertices, be traversed clockwise or counter-clockwise, and include or omit the collinear points lying in the middle of an edge. Without a fixed convention, an `==` would flag two correct results as different. That is why we fix: extreme points only, no duplicates, counter-clockwise order, starting at the lowest point (smallest `y`; ties broken by smallest `x`). With that, normalizing CGAL's output (with a `std::rotate`) is enough to compare thousands of cases without checking them by hand.

---

## Part 1 · The orientation test and Graham scan

### Implementation
- **`orient(a, b, c)`** computes `(bx − ax)(cy − ay) − (by − ay)(cx − ax)` using only `long long`. It returns the full value, not just its sign, because its magnitude (twice the area of triangle abc) is used by the rotating calipers in Part 5.
- **`graham(P)`**:
  1. `sort` + `unique` + `erase` with `operator<` to remove duplicates. If fewer than 3 points remain, they are returned as they are (since they are sorted, the lowest point comes first).
  2. Because of that ordering, `P[0]` is already the pivot p0 (lowest point; ties broken by leftmost).
  3. The remaining points are sorted by angle around p0 with a comparator built on `orient`: `a` comes before `b` if `orient(p0, a, b) > 0`; if they are collinear, the closer one comes first (`dist2`).
  4. Scan with a stack `H`: before pushing `p`, pop the top while `orient(H[-2], H[-1], p) ≤ 0`. The `≤` also removes collinear points lying in the middle of an edge.

All 4 test files (`example13`, `square_degenerate`, `collinear` and `single`) produce exactly the expected output.

### Question 1.1 · Can `orient` overflow?
With |x|, |y| ≤ 10⁹:

| Term | Maximum value |
|---|---|
| Each difference, e.g. `bx − ax` | 10⁹ − (−10⁹) = **2·10⁹** |
| Each product, e.g. `(bx − ax)(cy − ay)` | (2·10⁹)² = **4·10¹⁸** |
| The result (difference of two products) | 4·10¹⁸ + 4·10¹⁸ = **8·10¹⁸** |

A signed 64-bit `long long` goes up to 2⁶³ − 1 ≈ **9.22·10¹⁸**. Since 8·10¹⁸ < 9.22·10¹⁸, neither the intermediate products nor the result overflow (although the margin is small).

If coordinates could reach **10¹⁰**, differences would be up to 2·10¹⁰ and products up to 4·10²⁰, which **does** overflow 64 bits. The fix would be to compute the products in `__int128` (up to ~1.7·10³⁸) and compare only the sign. Using `double` would **not** work, because it loses exactness: it only represents integers exactly up to 2⁵³ ≈ 9·10¹⁵.

### Question 1.2 · Which points are popped in the class example?
Angular order around p0 = (312,56): p0, p1, p2, …, p12. In this example it matches the numbering used in class.

| When processing | Popped | Why |
|---|---|---|
| p3 (720,170) | p2 (636,104) | p1 → p2 → p3 turns right |
| p5 (576,194) | p4 (612,152) | |
| p6 (492,218) | p5 (576,194) | |
| p9 (402,224) | p8 (378,146), then p7 (420,200) | two pops in a row |
| p10 (420,308) | p9 (402,224), then p6 (492,218) | two pops in a row |
| p12 (192,146) | p11 (288,170) | |

There are **8 pops** in total, in the order **p2, p4, p5, p8, p7, p9, p6, p11**. The stack ends with **p0, p1, p3, p10, p12** = (312,56) (684,86) (720,170) (420,308) (192,146), which is the expected hull.

**Why is the scan O(n) even though the `while` loop can pop many points at once?** By amortized analysis: every point is pushed onto the stack **exactly once** and popped **at most once**. So the total number of pops over the whole scan is ≤ n, even if a single step pops several. There are n pushes and ≤ n pops, so the scan costs O(n). The overall O(n log n) cost of Graham scan is dominated by the sort.

---

## Part 2 · Jarvis march

### Implementation
Same preprocessing as Graham (remove duplicates and handle n < 3). Starting from `start = P[0]` (the lowest point), each iteration pushes `current` onto the hull and searches for the next vertex: it takes a candidate `next` and, for every point `r ≠ current`, computes `orient(current, next, r)`:
- if it is **< 0** (r lies to the right of current → next), `r` becomes the new candidate;
- if it is **= 0** (collinear) and `r` is **farther** from `current`, `r` is also taken, so points in the middle of an edge are not reported.

The loop (`do … while`) stops when it returns to `start`. Every call to `orient` is counted and the total is stored in `*tests`. The outputs of Jarvis and Graham are identical on all 4 test files.

### Question 2
On `example13.txt` (n = 13, h = 5) our Jarvis march performs **60 orientation tests**:

```
Get-Content data\example13.txt | .\lab.exe jarvis   →   orientation tests = 60
```

Each of the **h** iterations (one per hull vertex) compares the candidate against the other **n − 1** points, so:

**tests = h · (n − 1) = 5 · 12 = 60**

**Worst case:** when **every point is extreme** (h = n), for example points on a parabola `y = x²` or on a circle. The algorithm then performs n · (n − 1) tests, i.e. **Θ(n²)**. For n = 60 000 (the benchmark parabola) that is ≈ 3.6·10⁹ orientation tests.

---

## Part 3 · Stress testing against CGAL

### Implementation
- **`cgal_hull`**: after `CGAL::convex_hull_2`, the result is rotated so that it starts at the lowest point:
  `std::rotate(H.begin(), std::min_element(H.begin(), H.end()), H.end());`
  (`min_element` uses our `operator<`, which defines "lowest, then leftmost").
- **Hard cases:** besides the random points, depending on the seed we add
  - `seed % 3 == 1` → **duplicates**: up to 5 existing points are inserted again;
  - `seed % 3 == 2` → **collinear points on an edge**: the points (−1000 + 500k, −10⁹) for k = 0..4, which lie on the bottom edge. Only the two endpoints should appear in the hull.

Result on 3000 inputs (1000 seeds × 3 distributions):

```
OK: 0 mismatches
```

### Question 3
**Bugs found:** none. Graham and Jarvis already handled both degenerate cases by design: duplicates are removed with `sort` + `unique` before starting, and collinear points are excluded by the `≤ 0` condition in Graham's stack and by the "keep the farthest point" rule in Jarvis.

**Why is it safe to convert the coordinates to `double` before calling CGAL?** A `double` has a 53-bit mantissa, so it represents **exactly** every integer with |k| ≤ 2⁵³ ≈ 9.007·10¹⁵. Our coordinates satisfy |x|, |y| ≤ 10⁹ < 2⁵³, so the conversion loses no information and CGAL receives exactly the same points.

**Which kernel do we use, and what does it mean?** `CGAL::Exact_predicates_inexact_constructions_kernel` (EPICK):
- **Exact predicates:** questions that return a discrete decision (such as the left/right/collinear orientation) always give the correct answer. CGAL first evaluates them with interval arithmetic on `double` and, only if the result is ambiguous, recomputes them with exact arithmetic (GMP/MPFR).
- **Inexact constructions:** operations that create new geometric objects (e.g. the intersection point of two segments) are computed in `double` and may have rounding error.

The convex hull only uses **predicates** (it never constructs new points; it only selects points from the input), so with EPICK the result is exact and also fast.
