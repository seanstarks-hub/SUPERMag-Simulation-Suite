#!/usr/bin/env python3
"""Run Python-level benchmarks for SUPERMag."""

import time
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from supermag.proximity import pair_amplitude, critical_temperature


def bench_pair_amplitude():
    n_iter = 1000
    start = time.perf_counter()
    for _ in range(n_iter):
        pair_amplitude(d_F=50.0, xi_F=2.0, n_points=10000)
    elapsed = time.perf_counter() - start
    print(f"pair_amplitude: {n_iter} iters, {elapsed:.3f}s total, {elapsed/n_iter*1000:.3f}ms/call")


def bench_critical_temp():
    d_F = np.linspace(0.5, 50.0, 1000)
    n_iter = 100
    start = time.perf_counter()
    for _ in range(n_iter):
        critical_temperature(Tc0=9.2, d_S=50.0, d_F_array=d_F,
                             E_ex=256.0, xi_S=38.0, xi_F=0.7)
    elapsed = time.perf_counter() - start
    print(f"critical_temp: {n_iter} iters x {len(d_F)} pts, {elapsed:.3f}s total, {elapsed/n_iter*1000:.3f}ms/call")


if __name__ == "__main__":
    print("SUPERMag Python Benchmarks")
    print("-" * 40)
    bench_pair_amplitude()
    bench_critical_temp()
