#!/usr/bin/env python3
"""Sycamore-like Random Circuit Sampling by sparse Givens rotations."""

import argparse
import cmath
import math
import random
import sys
import time
from collections import namedtuple

import numpy as np
import psutil
from scipy.stats import unitary_group


def random_su2():
    """Sample a random SU(2) matrix."""
    U = unitary_group.rvs(2)
    det = np.linalg.det(U)
    return U / np.exp(1j * np.angle(det) / 2)


def random_su4():
    """Sample a random SU(4) matrix."""
    U = unitary_group.rvs(4)
    det = np.linalg.det(U)
    return U / np.exp(1j * np.angle(det) / 4)


BirthRecord = namedtuple(
    'BirthRecord',
    ['step', 'parent', 'child', 'i', 'j', 'theta', 'phi', 'dagger', 'conj_flag'],
)


def make_local_index(global_index, q1, q2):
    """Map a global basis index to a local 2-qubit index in {0,1,2,3}."""
    b1 = (global_index >> q1) & 1
    b2 = (global_index >> q2) & 1
    return b1 | (b2 << 1)


def make_global_indices(base_index, q1, q2):
    """Return the four global indices for the block of qubits q1,q2 with other bits fixed."""
    return [
        base_index,
        base_index | (1 << q1),
        base_index | (1 << q2),
        base_index | (1 << q1) | (1 << q2),
    ]


def complex_givens_to_zero(x, y, eps=1e-16):
    """Compute Givens theta, phi to zero y while acting on [x, y]^T."""
    if abs(y) < eps:
        return 0.0, 0.0
    if abs(x) < eps:
        return math.pi / 2, 0.0

    theta = math.atan2(abs(y), abs(x))
    phi = math.atan2(y.imag, y.real) - math.atan2(x.imag, x.real) + math.pi
    return theta, phi


def _record_birth(births, birth_map, step, parent, child, i, j, theta, phi, dagger, conj_flag):
    if child in birth_map:
        return
    record = BirthRecord(step, parent, child, i, j, theta, phi, dagger, conj_flag)
    births.append(record)
    birth_map[child] = record


def _apply_rotation_pair(
    amplitudes,
    i,
    j,
    theta,
    phi,
    dagger=False,
    step=None,
    births=None,
    birth_map=None,
    i_global=None,
    j_global=None,
    eps=1e-16,
):
    c = math.cos(theta)
    s = math.sin(theta)
    e = cmath.exp(-1j * phi)

    a_i = amplitudes.get(i, 0.0)
    a_j = amplitudes.get(j, 0.0)
    if dagger:
        a_i_new = c * a_i + e * s * a_j
        a_j_new = -e.conjugate() * s * a_i + c * a_j
    else:
        a_i_new = c * a_i - e * s * a_j
        a_j_new = e.conjugate() * s * a_i + c * a_j

    was_i = i in amplitudes
    was_j = j in amplitudes

    if abs(a_i_new) < eps:
        amplitudes.pop(i, None)
    else:
        amplitudes[i] = a_i_new
    if abs(a_j_new) < eps:
        amplitudes.pop(j, None)
    else:
        amplitudes[j] = a_j_new

    if births is None or birth_map is None or step is None:
        return

    if not was_i and abs(a_i_new) >= eps:
        parent = j
        child = i
        parent_global = j_global if j_global is not None else j
        child_global = i_global if i_global is not None else i
        conj_flag = not dagger
        _record_birth(
            births,
            birth_map,
            step,
            parent_global,
            child_global,
            i_global if i_global is not None else i,
            j_global if j_global is not None else j,
            theta,
            phi,
            dagger,
            conj_flag,
        )
    if not was_j and abs(a_j_new) >= eps:
        parent = i
        child = j
        parent_global = i_global if i_global is not None else i
        child_global = j_global if j_global is not None else j
        conj_flag = dagger
        _record_birth(
            births,
            birth_map,
            step,
            parent_global,
            child_global,
            i_global if i_global is not None else i,
            j_global if j_global is not None else j,
            theta,
            phi,
            dagger,
            conj_flag,
        )


def apply_rotation_pair(amplitudes, i, j, theta, phi, step=None, births=None, birth_map=None, i_global=None, j_global=None):
    """Apply a local Givens pairwise rotation to two complex amplitudes."""
    _apply_rotation_pair(
        amplitudes,
        i,
        j,
        theta,
        phi,
        dagger=False,
        step=step,
        births=births,
        birth_map=birth_map,
        i_global=i_global,
        j_global=j_global,
    )


def apply_rotation_pair_dagger(amplitudes, i, j, theta, phi, step=None, births=None, birth_map=None, i_global=None, j_global=None):
    """Apply the adjoint of a local Givens pairwise rotation."""
    _apply_rotation_pair(
        amplitudes,
        i,
        j,
        theta,
        phi,
        dagger=True,
        step=step,
        births=births,
        birth_map=birth_map,
        i_global=i_global,
        j_global=j_global,
    )


def compute_birth_depths(births, birth_map):
    depth_cache = {}

    def depth(index):
        if index not in birth_map:
            return 0
        if index in depth_cache:
            return depth_cache[index]
        parent = birth_map[index].parent
        result = 1 + depth(parent)
        depth_cache[index] = result
        return result

    return max((depth(child) for child in birth_map), default=0)


def estimate_birth_memory_bytes(births):
    return sys.getsizeof(births) + sum(sys.getsizeof(rec) for rec in births)


def recover_zero_structure(births):
    return [
        {
            'step': rec.step,
            'parent': rec.parent,
            'child': rec.child,
            'i': rec.i,
            'j': rec.j,
            'theta': rec.theta,
            'phi': rec.phi,
            'dagger': rec.dagger,
        }
        for rec in births
    ]


def reconstruct_parent_amplitude_at_birth(child_amplitude, record):
    """Reconstruct the parent amplitude immediately before a birth event.

    This assumes `child_amplitude` is the amplitude of the child immediately
    after the Givens step that created it.
    """
    s = math.sin(record.theta)
    if abs(s) < 1e-16:
        raise ValueError('Cannot reconstruct parent amplitude with s=0.')
    if record.dagger:
        if record.child == record.i:
            return cmath.exp(1j * record.phi) * child_amplitude / s
        return -cmath.exp(-1j * record.phi) * child_amplitude / s
    if record.child == record.j:
        return cmath.exp(1j * record.phi) * child_amplitude / s
    return -cmath.exp(-1j * record.phi) * child_amplitude / s


def reconstruct_amplitude(index, step, births, birth_map, final_state=None):
    """Estimate the amplitude of a basis index at the given step using the birth tree.

    This function can determine that an index was zero before its birth step.
    Full amplitude recovery after the birth step requires the full evolution
    on that index and is not available from birth records alone.
    """
    if step <= 0:
        return 1.0 + 0j if index == 0 else 0.0
    record = birth_map.get(index)
    if record is None:
        if index == 0:
            return 1.0 + 0j
        if final_state is not None:
            return final_state.get(index, 0.0)
        return None
    if step < record.step:
        return 0.0
    if step == record.step:
        return 0.0
    return None


def left_mult_givens_rows(A, i, j, theta, phi):
    """Apply left-multiplication by a 2-level Givens on rows i and j."""
    c = math.cos(theta)
    s = math.sin(theta)
    e = cmath.exp(-1j * phi)

    row_i = A[i, :].copy()
    row_j = A[j, :].copy()
    A[i, :] = c * row_i - s * e * row_j
    A[j, :] = s * e.conjugate() * row_i + c * row_j


def right_mult_givens_columns(A, i, j, theta, phi):
    """Apply right-multiplication by a 2-level Givens on columns i and j."""
    c = math.cos(theta)
    s = math.sin(theta)
    e = cmath.exp(-1j * phi)

    col_i = A[:, i].copy()
    col_j = A[:, j].copy()
    A[:, i] = c * col_i - s * e * col_j
    A[:, j] = s * e.conjugate() * col_i + c * col_j


def decompose_su4_to_givens(U, q1, q2, n_qubits):
    """Decompose an SU(4) two-qubit unitary into local 2-level Givens rotations.

    Returns a tuple (givens_ops, diagonal_phases), where givens_ops is a list of
    (local_i, local_j, theta, phi) and diagonal_phases is a length-4 complex array.
    """
    assert U.shape == (4, 4)
    assert U.dtype == np.complex128 or U.dtype == np.complex64

    A = U.copy().astype(np.complex128)
    givens_ops = []

    # Zero lower-triangular entries with left-multiplying Givens.
    # The order below zeros (3,0), (2,0), (1,0), (3,1), (2,1), (3,2).
    step_order = [
        ((2, 3), 0),
        ((1, 2), 0),
        ((0, 1), 0),
        ((2, 3), 1),
        ((1, 2), 1),
        ((2, 3), 2),
    ]

    for (i, j), col in step_order:
        x = A[i, col]
        y = A[j, col]
        theta, phi = complex_givens_to_zero(x, y)
        if abs(theta) > 1e-15 or abs(phi) > 1e-15:
            left_mult_givens_rows(A, i, j, theta, phi)
        givens_ops.append((i, j, theta, phi))

    diag_phases = np.array([1.0 + 0j] * 4, dtype=np.complex128)
    for k in range(4):
        if abs(A[k, k]) < 1e-16:
            diag_phases[k] = 1.0 + 0j
        else:
            diag_phases[k] = A[k, k] / abs(A[k, k])

    # Remove the diagonal phase from R to make the remaining multiplier unitary.
    # The decomposition U = (G_1^† ... G_6^†) diag(diag_phases)
    return givens_ops, diag_phases


def apply_local_givens_to_state(
    state,
    q1,
    q2,
    givens_ops,
    diag_phases=None,
    track=False,
    births=None,
    birth_map=None,
    step_counter=1,
):
    """Apply local 4D Givens operations to a global sparse state."""
    if q1 > q2:
        q1, q2 = q2, q1

    other_mask = ~((1 << q1) | (1 << q2))
    block_map = {}

    # Group amplitudes by the fixed bits outside q1,q2.
    for index, amplitude in state.items():
        base = index & other_mask
        local_index = (((index >> q1) & 1) << 1) | ((index >> q2) & 1)
        block = block_map.setdefault(base, {})
        block[local_index] = amplitude

    if diag_phases is not None:
        for block in block_map.values():
            for local_index, phase in enumerate(diag_phases):
                if local_index in block:
                    block[local_index] *= phase

    history = []
    for local_i, local_j, theta, phi in reversed(givens_ops):
        for base, block in block_map.items():
            i_global = base | (((local_i >> 1) & 1) << q1) | ((local_i & 1) << q2)
            j_global = base | (((local_j >> 1) & 1) << q1) | ((local_j & 1) << q2)
            apply_rotation_pair_dagger(
                block,
                local_i,
                local_j,
                theta,
                phi,
                step=step_counter,
                births=births,
                birth_map=birth_map,
                i_global=i_global,
                j_global=j_global,
            )
        if track:
            history.append(sum(len(block) for block in block_map.values()))
        step_counter += 1

    state.clear()
    for base, block in block_map.items():
        for local_index, amplitude in block.items():
            if abs(amplitude) < 1e-16:
                continue
            q1_bit = (local_index >> 1) & 1
            q2_bit = local_index & 1
            global_index = base | (q1_bit << q1) | (q2_bit << q2)
            state[global_index] = amplitude

    return (history if track else None), step_counter


def decompose_su2_to_givens(U, qubit, n_qubits):
    """Decompose an SU(2) single-qubit unitary into a local 2-level Givens rotation.

    For SU(2), a single complex Givens plus local phases is sufficient.
    """
    assert U.shape == (2, 2)
    assert abs(np.linalg.det(U) - 1) < 1e-8

    A = U.copy().astype(np.complex128)
    x = A[0, 0]
    y = A[1, 0]
    theta, phi = complex_givens_to_zero(x, y)
    left_mult_givens_rows(A, 0, 1, theta, phi)

    diag_phases = np.array([1.0 + 0j, 1.0 + 0j], dtype=np.complex128)
    for k in range(2):
        if abs(A[k, k]) > 1e-16:
            diag_phases[k] = A[k, k] / abs(A[k, k])

    return [(0, 1, theta, phi)], diag_phases


def apply_single_qubit_gate(
    state,
    qubit,
    givens_ops,
    diag_phases=None,
    track=False,
    births=None,
    birth_map=None,
    step_counter=1,
):
    """Apply a single-qubit gate to the global sparse state using local rotations."""
    other_mask = ~(1 << qubit)
    block_map = {}

    for index, amplitude in state.items():
        base = index & other_mask
        local_index = (index >> qubit) & 1
        block = block_map.setdefault(base, {})
        block[local_index] = amplitude

    if diag_phases is not None:
        for block in block_map.values():
            for local_index, phase in enumerate(diag_phases):
                if local_index in block:
                    block[local_index] *= phase

    history = []
    for local_i, local_j, theta, phi in reversed(givens_ops):
        for base, block in block_map.items():
            i_global = base | (local_i << qubit)
            j_global = base | (local_j << qubit)
            apply_rotation_pair_dagger(
                block,
                local_i,
                local_j,
                theta,
                phi,
                step=step_counter,
                births=births,
                birth_map=birth_map,
                i_global=i_global,
                j_global=j_global,
            )
        if track:
            history.append(sum(len(block) for block in block_map.values()))
        step_counter += 1

    state.clear()
    for base, block in block_map.items():
        for local_index, amplitude in block.items():
            if abs(amplitude) < 1e-16:
                continue
            global_index = base | (local_index << qubit)
            state[global_index] = amplitude

    return (history if track else None), step_counter


def build_sycamore_like_circuit(n_qubits, depth, seed=None):
    """Generate a Sycamore-like circuit with alternating 1q and 2q layers."""
    random.seed(seed)
    circuit = []

    for layer in range(depth):
        # Single-qubit random SU(2) on all qubits.
        for q in range(n_qubits):
            circuit.append(('one_qubit', q, random_su2()))

        # Two-qubit layer with staggered nearest-neighbour couplings.
        if layer % 2 == 0:
            pairs = [(q, q + 1) for q in range(0, n_qubits - 1, 2)]
        else:
            pairs = [(q, q + 1) for q in range(1, n_qubits - 1, 2)]

        for q1, q2 in pairs:
            circuit.append(('two_qubit', q1, q2, random_su4()))

    return circuit


def sample_state(state, n_samples=10):
    """Sample bitstrings from the sparse quantum state."""
    probs = [abs(amplitude) ** 2 for amplitude in state.values()]
    total = sum(probs)
    if total == 0:
        return []
    normalized = [p / total for p in probs]
    indices = list(state.keys())

    samples = []
    for _ in range(n_samples):
        r = random.random()
        cum = 0.0
        for idx, p in zip(indices, normalized):
            cum += p
            if r <= cum:
                samples.append(idx)
                break
    return samples


def format_bitstring(index, n_qubits):
    return format(index, '0{}b'.format(n_qubits))


def generate_random_givens_circuit(n_qubits, n_rotations, seed=None):
    """Generate a random sequence of global Givens rotations."""
    random.seed(seed)
    circuit = []

    for _ in range(n_rotations):
        i = random.getrandbits(n_qubits)
        j = random.getrandbits(n_qubits)
        if j >= i:
            j += 1
        theta = random.uniform(0, math.pi / 2)
        phi = random.uniform(0, 2 * math.pi)
        circuit.append((i, j, theta, phi))

    return circuit


def apply_global_givens(state, i, j, theta, phi, step=None, births=None, birth_map=None):
    """Apply a global Givens rotation directly to the sparse state."""
    apply_rotation_pair(
        state,
        i,
        j,
        theta,
        phi,
        step=step,
        births=births,
        birth_map=birth_map,
        i_global=i,
        j_global=j,
    )


def verify_support_history(history):
    return all(support <= 1 + k for k, support in enumerate(history, start=1))


def run_givens(n_qubits, n_rotations, seed=42, sample_count=5, track_births=False):
    circuit = generate_random_givens_circuit(n_qubits, n_rotations, seed)
    state = {0: 1.0}
    support_history = []
    births = [] if track_births else None
    birth_map = {} if track_births else None
    step_counter = 1
    process = psutil.Process()
    peak_rss = process.memory_info().rss

    start = time.perf_counter()
    for i, j, theta, phi in circuit:
        apply_global_givens(
            state,
            i,
            j,
            theta,
            phi,
            step=step_counter,
            births=births,
            birth_map=birth_map,
        )
        support_history.append(len(state))
        step_counter += 1
        rss = process.memory_info().rss
        if rss > peak_rss:
            peak_rss = rss

    elapsed = time.perf_counter() - start
    samples = sample_state(state, sample_count)
    bitstrings = [format_bitstring(idx, n_qubits) for idx in samples]
    birth_depth = compute_birth_depths(births, birth_map) if births is not None else None
    birth_bytes = estimate_birth_memory_bytes(births) if births is not None else None
    example_birth = recover_zero_structure(births[:1])[0] if births else None
    return {
        'mode': 'givens',
        'n_qubits': n_qubits,
        'n_rotations': n_rotations,
        'total_givens': len(circuit),
        'final_support': len(state),
        'birth_count': len(births) if births is not None else None,
        'birth_depth': birth_depth,
        'birth_bytes': birth_bytes,
        'birth_count_ok': len(births) == max(0, len(state) - 1) if births is not None else None,
        'example_birth': example_birth,
        'max_support': max(support_history, default=1),
        'support_bound_ok': verify_support_history(support_history),
        'time_s': elapsed,
        'peak_rss_mb': peak_rss / (1024**2),
        'samples': bitstrings,
    }


def run_rcs(n_qubits, depth, seed=42, sample_count=5, track_births=False):
    circuit = build_sycamore_like_circuit(n_qubits, depth, seed)
    state = {0: 1.0}
    total_givens = 0
    total_one_qubit = 0
    total_two_qubit = 0
    support_history = []
    births = [] if track_births else None
    birth_map = {} if track_births else None
    step_counter = 1
    process = psutil.Process()
    peak_rss = process.memory_info().rss

    start = time.perf_counter()
    for gate in circuit:
        if gate[0] == 'one_qubit':
            _, q, U = gate
            givens_ops, diag_phases = decompose_su2_to_givens(U, q, n_qubits)
            history, step_counter = apply_single_qubit_gate(
                state,
                q,
                givens_ops,
                diag_phases,
                track=True,
                births=births,
                birth_map=birth_map,
                step_counter=step_counter,
            )
            total_one_qubit += 1
        else:
            _, q1, q2, U = gate
            givens_ops, diag_phases = decompose_su4_to_givens(U, q1, q2, n_qubits)
            history, step_counter = apply_local_givens_to_state(
                state,
                q1,
                q2,
                givens_ops,
                diag_phases,
                track=True,
                births=births,
                birth_map=birth_map,
                step_counter=step_counter,
            )
            total_two_qubit += 1

        for support in history:
            total_givens += 1
            support_history.append(support)

        rss = process.memory_info().rss
        if rss > peak_rss:
            peak_rss = rss

    elapsed = time.perf_counter() - start
    samples = sample_state(state, sample_count)
    bitstrings = [format_bitstring(idx, n_qubits) for idx in samples]
    birth_depth = compute_birth_depths(births, birth_map) if births is not None else None
    birth_bytes = estimate_birth_memory_bytes(births) if births is not None else None
    example_birth = recover_zero_structure(births[:1])[0] if births else None
    return {
        'mode': 'decompose',
        'n_qubits': n_qubits,
        'depth': depth,
        'total_gates': len(circuit),
        'one_qubit_gates': total_one_qubit,
        'two_qubit_gates': total_two_qubit,
        'total_givens': total_givens,
        'final_support': len(state),
        'birth_count': len(births) if births is not None else None,
        'birth_depth': birth_depth,
        'birth_bytes': birth_bytes,
        'birth_count_ok': len(births) == max(0, len(state) - 1) if births is not None else None,
        'example_birth': example_birth,
        'max_support': len(state),
        'support_bound_ok': None,
        'time_s': elapsed,
        'peak_rss_mb': peak_rss / (1024**2),
        'samples': bitstrings,
    }


def parse_args():
    parser = argparse.ArgumentParser(description='Simulate Sycamore-like RCS using sparse Givens rotations.')
    parser.add_argument('--mode', choices=['decompose', 'givens'], default='givens', help='Simulation mode: givens = direct 2-level Givens theorem demo, decompose = SU(4) decomposition into local Givens')
    parser.add_argument('--n-qubits', type=int, default=10, help='Number of qubits')
    parser.add_argument('--depth', type=int, default=4, help='Number of layers for decompose mode')
    parser.add_argument('--n-rotations', type=int, default=100, help='Number of Givens rotations for givens mode')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--samples', type=int, default=5, help='Number of sampled bitstrings')
    parser.add_argument('--track-births', action='store_true', help='Collect birth tree records during simulation')
    return parser.parse_args()


def print_result(result):
    print('mode=', result['mode'])
    print('n_qubits=', result['n_qubits'])
    if result['mode'] == 'decompose':
        print('depth=', result['depth'])
        print('total_gates=', result['total_gates'])
        print('one_qubit_gates=', result['one_qubit_gates'])
        print('two_qubit_gates=', result['two_qubit_gates'])
    else:
        print('n_rotations=', result['n_rotations'])
    print('total_givens=', result['total_givens'])
    print('final_support=', result['final_support'])
    print('birth_count=', 'N/A' if result.get('birth_count') is None else result['birth_count'])
    print('birth_depth=', 'N/A' if result.get('birth_depth') is None else result['birth_depth'])
    birth_bytes = result.get('birth_bytes')
    print('birth_bytes=', 'N/A' if birth_bytes is None else f"{birth_bytes:,}")
    print('birth_count_ok=', 'N/A' if result.get('birth_count_ok') is None else result['birth_count_ok'])
    if 'example_birth' in result and result['example_birth'] is not None:
        print('example_birth=', result['example_birth'])
    print('max_support=', result['max_support'])
    print('support_bound_ok=', 'N/A' if result['support_bound_ok'] is None else result['support_bound_ok'])
    print('time_s=', f"{result['time_s']:.6f}")
    print('peak_rss_mb=', f"{result['peak_rss_mb']:.3f}")
    print('samples=')
    for bitstring in result['samples']:
        print('  ', bitstring)


def main():
    args = parse_args()
    if args.mode == 'decompose':
        result = run_rcs(args.n_qubits, args.depth, args.seed, args.samples, track_births=args.track_births)
    else:
        result = run_givens(args.n_qubits, args.n_rotations, args.seed, args.samples, track_births=args.track_births)
    print_result(result)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
