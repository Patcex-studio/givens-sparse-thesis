#!/usr/bin/env python3
"""
Verification script for Givens Sparse Thesis:
Demonstrates that the support bound |supp(ψ_k)| ≤ 1 + k holds.
"""

import cmath
import random
import math

def givens_rotation(theta, phi, a_i, a_j):
    """
    Apply Givens rotation to amplitudes a_i and a_j.
    
    G(θ,φ) = [[c, -e^{-iφ}s], [e^{iφ}s, c]]
    where c = cosθ, s = sinθ
    """
    c = cmath.cos(theta)
    s = cmath.sin(theta)
    e = cmath.exp(-1j * phi)
    
    # Forward update
    a_i_new = c * a_i - e * s * a_j
    a_j_new = e.conjugate() * s * a_i + c * a_j
    
    return a_i_new, a_j_new

def remove_zeros(state, eps=1e-12):
    """Remove keys with amplitudes close to zero."""
    to_del = [k for k, v in state.items() if abs(v) < eps]
    for k in to_del:
        del state[k]

def apply_rotation(state, i, j, theta, phi):
    """Apply Givens rotation to the quantum state."""
    if i in state and j in state:
        a_i = state[i]
        a_j = state[j]
        state[i], state[j] = givens_rotation(theta, phi, a_i, a_j)
    elif i in state:
        a_i = state[i]
        a_j = 0.0
        state[i], state[j] = givens_rotation(theta, phi, a_i, a_j)
    elif j in state:
        a_j = state[j]
        a_i = 0.0
        state[i], state[j] = givens_rotation(theta, phi, a_i, a_j)
    # If neither, both remain zero
    remove_zeros(state)

def verify_support_bound(rotations, n_qubits=3):
    """
    Verify the sparse support bound: |supp(ψ_k)| ≤ 1 + k
    
    Args:
        rotations: list of (i, j, theta, phi) tuples
        n_qubits: number of qubits (determines basis size)
    
    Returns:
        True if bound holds, False otherwise
    """
    basis_size = 2 ** n_qubits
    state = {0: 1.0}  # Start with |0...0⟩
    
    # Initial state (before any rotations)
    initial_support = len(state)
    
    for k, (i, j, theta, phi) in enumerate(rotations):
        # Before applying rotation k+1, we have done k rotations
        # The bound after k rotations is: |supp| ≤ 1 + k
        support_size_before = len(state)
        max_allowed_before = 1 + k
        
        if support_size_before > max_allowed_before:
            print(f"VIOLATION before step {k+1}: support = {support_size_before}, max allowed = {max_allowed_before}")
            return False
        
        apply_rotation(state, i, j, theta, phi)
        support_size_after = len(state)
        
        # After applying rotation, we've done k+1 rotations
        # The new bound is: |supp| ≤ 1 + (k+1)
        max_allowed_after = 1 + (k + 1)
        
        if support_size_after > max_allowed_after:
            print(f"VIOLATION after step {k+1}: support = {support_size_after}, max allowed = {max_allowed_after}")
            return False
    
    return True

def main():
    print("=" * 60)
    print("Givens Sparse Thesis Verification")
    print("=" * 60)
    
    # Example 1: Simple circuit with 3 rotations on 2 qubits
    print("\n[Example 1] 3 rotations on 2 qubits")
    print("-" * 40)
    
    rotations1 = [
        (0, 1, math.pi/4, 0),       # G(0,1,π/4,0)
        (0, 1, math.pi/4, math.pi/2),  # G(0,1,π/4,π/2)
        (0, 1, math.pi/4, 0),       # G(0,1,π/4,0)
    ]
    
    print("Rotations applied:")
    for k, (i, j, t, p) in enumerate(rotations1):
        print(f"  Step {k+1}: G({i},{j}, {t:.4f}, {p:.4f})")
    
    result1 = verify_support_bound(rotations1, n_qubits=2)
    print(f"\nSupport bound verified: {result1}")
    
    # Show final state
    basis_size = 2 ** 2
    state = {0: 1.0}
    for (i, j, t, p) in rotations1:
        apply_rotation(state, i, j, t, p)
    print(f"Final state support: {len(state)} (max allowed: {1 + len(rotations1)})")
    
    # Example 2: Random rotations
    print("\n[Example 2] 10 random rotations on 3 qubits")
    print("-" * 40)
    
    random.seed(42)
    rotations2 = []
    for k in range(10):
        i = random.randint(0, 7)
        j = random.randint(0, 7)
        if i == j:
            j = (j + 1) % 8
        theta = random.uniform(0, 2 * math.pi / 4)
        phi = random.uniform(0, 2 * math.pi)
        rotations2.append((i, j, theta, phi))
    
    print(f"Applied {len(rotations2)} random rotations")
    result2 = verify_support_bound(rotations2, n_qubits=3)
    print(f"Support bound verified: {result2}")
    
    # Example 3: Stress test with many rotations
    print("\n[Example 3] Stress test: 100 rotations on 4 qubits")
    print("-" * 40)
    
    random.seed(123)
    rotations3 = []
    for k in range(100):
        i = random.randint(0, 15)
        j = random.randint(0, 15)
        if i == j:
            j = (j + 1) % 16
        theta = random.uniform(0, 2 * math.pi / 4)
        phi = random.uniform(0, 2 * math.pi)
        rotations3.append((i, j, theta, phi))
    
    print(f"Applied {len(rotations3)} random rotations")
    result3 = verify_support_bound(rotations3, n_qubits=4)
    print(f"Support bound verified: {result3}")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Example 1 (3 rotations): {'PASS' if result1 else 'FAIL'}")
    print(f"Example 2 (10 rotations): {'PASS' if result2 else 'FAIL'}")
    print(f"Example 3 (100 rotations): {'PASS' if result3 else 'FAIL'}")
    
    if result1 and result2 and result3:
        print("\n✓ All tests passed! The sparse support bound holds.")
    else:
        print("\n✗ Some tests failed!")
    
    return result1 and result2 and result3

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)