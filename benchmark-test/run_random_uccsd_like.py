#!/usr/bin/env python3
"""Run random UCCSD-like Givens rotations and measure sparse support."""

import argparse
import math
import random
import time

import psutil


def remove_zeros(state, eps=1e-12):
    for key in [k for k, v in state.items() if abs(v) < eps]:
        del state[key]


def apply_rotation(state, i, j, theta, phi):
    c = math.cos(theta)
    s = math.sin(theta)
    e = complex(math.cos(-phi), math.sin(-phi))

    if i in state and j in state:
        a_i = state[i]
        a_j = state[j]
        state[i] = c * a_i - e * s * a_j
        state[j] = e.conjugate() * s * a_i + c * a_j
    elif i in state:
        a_i = state[i]
        state[i] = c * a_i
        state[j] = e.conjugate() * s * a_i
    elif j in state:
        a_j = state[j]
        state[i] = -e * s * a_j
        state[j] = c * a_j
    remove_zeros(state)


def generate_uccsd_like_rotation(n_qubits):
    n_bits = random.choice([1, 2])
    i = random.getrandbits(n_qubits)

    if n_bits == 1:
        bit = random.randrange(n_qubits)
        j = i ^ (1 << bit)
    else:
        b1, b2 = random.sample(range(n_qubits), 2)
        j = i ^ ((1 << b1) | (1 << b2))

    theta = random.uniform(0, 2 * math.pi)
    phi = random.uniform(0, 2 * math.pi)
    return i, j, theta, phi


def run_experiment(n_qubits, m, seed=42):
    random.seed(seed)
    state = {0: 1.0}
    process = psutil.Process()
    peak_rss = process.memory_info().rss
    max_support = len(state)
    start = time.perf_counter()

    for _ in range(m):
        i, j, theta, phi = generate_uccsd_like_rotation(n_qubits)
        apply_rotation(state, i, j, theta, phi)
        support = len(state)
        if support > max_support:
            max_support = support
        rss = process.memory_info().rss
        if rss > peak_rss:
            peak_rss = rss

    end = time.perf_counter()
    return len(state), max_support, end - start, peak_rss


def parse_args():
    parser = argparse.ArgumentParser(description='Run random UCCSD-like sparse support experiment.')
    parser.add_argument('--n-qubits', type=int, default=148, help='Number of qubits')
    parser.add_argument('--m', type=int, default=1000, help='Number of rotations')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    return parser.parse_args()


def main():
    args = parse_args()
    final_support, max_support, elapsed, peak_rss = run_experiment(args.n_qubits, args.m, args.seed)

    print('n_qubits=', args.n_qubits)
    print('M=', args.m)
    print('final_support=', final_support)
    print('max_support=', max_support)
    print('time_s=', f'{elapsed:.6f}')
    print('peak_rss_mb=', f'{peak_rss / (1024*1024):.3f}')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
