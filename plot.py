# Part 4: plots from results/results.csv (made with: ./bench > results/results.csv)
# Usage: python3 plot.py   ->  writes the PNGs into results/
import csv
import math
import os
import matplotlib.pyplot as plt

DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
DISTS = ["square", "disk", "parabola"]
ALGOS = ["graham", "jarvis", "cgal_default", "cgal_graham_andrew", "cgal_jarvis"]

with open(os.path.join(DIR, "results.csv")) as f:
    rows = list(csv.DictReader(f))

# (a) h vs n, one line per distribution. Every algorithm repeats the same h,
# so we only read the graham rows.
plt.figure()
for d in DISTS:
    pts = [(int(r["n"]), int(r["h"])) for r in rows if r["dist"] == d and r["algo"] == "graham"]
    ns = [p[0] for p in pts]
    hs = [p[1] for p in pts]
    plt.loglog(ns, hs, "o-", label=d)
    # on a log-log plot n^k is a straight line with slope k, so this estimates k
    slope = math.log(hs[-1] / hs[0]) / math.log(ns[-1] / ns[0])
    print(f"slope of h vs n for {d}: {slope:.3f}")
plt.xlabel("n (points)")
plt.ylabel("h (hull vertices)")
plt.title("Hull size h vs n")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.savefig(os.path.join(DIR, "h_vs_n.png"), dpi=150)

# (b) time vs n, one plot per distribution and one line per algorithm
for d in DISTS:
    plt.figure()
    for a in ALGOS:
        pts = [(int(r["n"]), float(r["ms"])) for r in rows if r["dist"] == d and r["algo"] == a]
        if pts:   # the Jarvis rows can be missing if bench skipped them
            plt.loglog([p[0] for p in pts], [p[1] for p in pts], "o-", label=a)
    plt.xlabel("n (points)")
    plt.ylabel("time (ms)")
    plt.title(f"Running time vs n ({d})")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.savefig(os.path.join(DIR, f"time_{d}.png"), dpi=150)

print("plots saved in", DIR)
