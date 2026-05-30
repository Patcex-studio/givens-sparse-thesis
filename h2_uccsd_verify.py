#!/usr/bin/env python3
"""H2 UCCSD verification using OpenFermion test data.

This script loads the built-in OpenFermion H2/STO-3G dataset, builds a
4-qubit UCCSD ansatz in the (2e, 2o) active space, optimizes the single
relevant amplitude, and checks support growth for a direct Givens rotation.
"""

import cmath
import math
import os
import time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from openfermion.chem import MolecularData
from openfermion.config import DATA_DIRECTORY
from openfermion.circuits import uccsd_singlet_generator
from openfermion.linalg.sparse_tools import get_sparse_operator, jw_hartree_fock_state
from openfermion.transforms.opconversions import get_fermion_operator, jordan_wigner, normal_ordered


def apply_rotation(state, i, j, theta, phi):
    """Apply a Givens rotation to sparse state amplitudes."""
    e = cmath.exp(-1j * phi)
    c = math.cos(theta)
    s = math.sin(theta)

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


def remove_zeros(state, eps=1e-12):
    for k in [k for k, v in state.items() if abs(v) < eps]:
        del state[k]


def extract_givens_angles(a_i, a_j, eps=1e-12):
    """Recover Givens rotation parameters from the final two-state amplitudes."""
    if abs(a_i.imag) > eps:
        raise ValueError('Expected real amplitude for the initial basis state i.')

    c = a_i.real
    if c > 1.0:
        c = 1.0
    elif c < -1.0:
        c = -1.0

    theta = math.acos(c)
    s = math.sin(theta)
    if abs(s) < eps:
        phi = 0.0
    else:
        phi = -math.atan2(a_j.imag, a_j.real)
        if phi < 0.0:
            phi += 2 * math.pi

    return theta, phi


def build_h2_hamiltonian():
    filename = os.path.join(DATA_DIRECTORY, 'H2_sto-3g_singlet_0.7414')
    molecule = MolecularData(
        [('H', (0.0, 0.0, 0.0)), ('H', (0.0, 0.0, 0.7414))],
        'sto-3g',
        1,
        filename=filename,
    )
    molecule.load()
    molecular_hamiltonian = molecule.get_molecular_hamiltonian()
    fermion_hamiltonian = normal_ordered(get_fermion_operator(molecular_hamiltonian))
    qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
    qubit_matrix = get_sparse_operator(qubit_hamiltonian).toarray()
    return molecule, qubit_matrix


def energy_for_theta(theta, hf, qubit_matrix):
    generator = uccsd_singlet_generator([0.0, theta], 4, 2)
    U = expm(get_sparse_operator(generator).toarray())
    psi = U.dot(hf)
    return np.vdot(psi, qubit_matrix.dot(psi)).real, psi


def main():
    t0 = time.perf_counter()
    molecule, qubit_matrix = build_h2_hamiltonian()
    hf = jw_hartree_fock_state(2, 4)

    result = minimize_scalar(
        lambda theta: energy_for_theta(theta, hf, qubit_matrix)[0],
        bounds=(-2.0, 2.0),
        method='bounded',
        options={'xatol': 1e-10},
    )

    theta_opt = result.x
    energy_opt, psi_opt = energy_for_theta(theta_opt, hf, qubit_matrix)

    # Support of the final state.
    final_support = int(np.count_nonzero(np.abs(psi_opt) > 1e-12))

    # This ansatz is a single two-level rotation between |1100> (12) and |0011> (3).
    i, j = 12, 3
    a_i = psi_opt[i]
    a_j = psi_opt[j]
    theta_givens, phi_givens = extract_givens_angles(a_i, a_j)

    # Reconstruct with a single Givens rotation for verification.
    state = {i: 1.0}
    apply_rotation(state, i, j, theta_givens, phi_givens)
    final_support_givens = len(state)
    support_path = max(1, final_support_givens)

    t1 = time.perf_counter()
    runtime = t1 - t0

    print('H2 UCCSD: n_qubits=4, M=1, final |supp|=%d, time=%.6f sec, energy=%.12f (совпадает? %s)' % (
        final_support,
        runtime,
        energy_opt,
        'да' if abs(energy_opt - molecule.fci_energy) < 1e-9 else 'нет',
    ))
    print('Hartree-Fock index:', np.nonzero(np.abs(hf) > 1e-12)[0])
    print('Ansatz nonzero indices:', np.nonzero(np.abs(psi_opt) > 1e-12)[0])
    print('Single Givens rotation: i=%d, j=%d, theta=%.12f, phi=%.12f' % (i, j, theta_givens, phi_givens))
    print('Givens final support:', final_support_givens)
    print('FCI energy:', molecule.fci_energy)

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
