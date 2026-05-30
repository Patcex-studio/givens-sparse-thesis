# Exact Classical Simulation of Givens-Circuit Quantum States

**Abstract**  
We present a linear-time, exact classical algorithm for simulating quantum circuits composed solely of complex Givens rotations. The key insight is that each Givens rotation affects exactly two basis states, leading to a linear bound on the number of non-zero amplitudes: after M rotations, at most 1+M amplitudes are non-zero. This yields O(M) time and memory complexity, avoiding the exponential blowup that typically plagues quantum circuit simulation.

## 1. Introduction

Quantum circuit simulation on classical computers faces an exponential barrier: a circuit with n qubits requires O(2^n) memory to represent the full state vector. This has been cited as evidence that quantum computers are necessary for simulating quantum systems. However, this barrier assumes worst-case circuits that explore the entire Hilbert space.

We identify a non-trivial class of circuits — those built from complex Givens rotations — that admit efficient classical simulation. A Givens rotation G(i,j,θ,φ) acts on exactly two computational basis states |i⟩ and |j⟩, performing a rotation in the two-dimensional subspace spanned by these states.

## 2. The Sparse Support Theorem

**Theorem 1 (Sparse Support Bound).** For any sequence of M complex Givens rotations applied to the initial state |0…0⟩, the number of non-zero amplitudes satisfies:

$$|\text{supp}(\psi_k)| \le 1 + k \quad \text{for all } k = 0,\dots,M$$

*Proof.* By induction. Base case k=0: ψ₀ = |0⟩ has exactly one non-zero amplitude. Inductive step: when applying G(i,j,θ,φ), at most one new basis state becomes non-zero. Thus the support grows by at most 1 per rotation. ∎

This linear bound implies that for circuits with O(n) rotations, the state can be represented in O(n) memory — polynomial in the number of qubits.

## 3. Linear-Time Algorithm

### 3.1 Data Structure

Our simulator uses a hash map (dictionary) storing only non-zero amplitudes:

$$\text{map}: \text{basis index} \to \text{complex amplitude}$$

### 3.2 Forward Pass

For each rotation G(i,j,θ,φ), we perform the following update:

```
initialize map[0] = 1.0
for each rotation (i, j, θ, φ):
    a_i = map.get(i, 0)
    a_j = map.get(j, 0)
    c = cosθ, s = sinθ, e = exp(-iφ)
    a_i' = c·a_i - e^{-iφ}·s·a_j
    a_j' = e^{iφ}·s·a_i + c·a_j
    map[i] = a_i'
    map[j] = a_j'
```

Each step requires O(1) dictionary operations. Total forward pass: O(M).

### 3.3 Expectation Values

For a Pauli string P with mask m and phase z, the expectation value is:

$$\langle\psi|P|\psi\rangle = \sum_{x\in\text{supp}} \psi_x^* \psi_{x\oplus m} (-1)^{x\cdot z} \cdot i^{|m\wedge z|}$$

The sum runs over O(|supp|) = O(M) terms. For L observables, total time: O(L·M).

## 4. Backpropagation for Gradient Computation

### 4.1 Rotation Matrix and Adjoint

The Givens rotation matrix is:

$$G(\theta,\phi) = \begin{bmatrix} c & -e^{-i\phi}s \\ e^{i\phi}s & c \end{bmatrix}$$

where c = cosθ, s = sinθ.

Its adjoint (inverse) is:

$$G^\dagger(\theta,\phi) = \begin{bmatrix} c & e^{-i\phi}s \\ -e^{i\phi}s & c \end{bmatrix}$$

### 4.2 Forward Update

For amplitudes a_i and a_j before rotation:

$$\begin{aligned}
a_i' &= c\cdot a_i - e^{-i\phi}s\cdot a_j \\
a_j' &= e^{i\phi}s\cdot a_i + c\cdot a_j
\end{aligned}$$

### 4.3 Backward (Reverse Pass)

Given post-rotation amplitudes a_i', a_j' and adjoint entries b_i, b_j, we first restore pre-rotation amplitudes:

$$\begin{aligned}
a_i &= c\cdot a_i' + e^{-i\phi}s\cdot a_j' \\
a_j &= -e^{i\phi}s\cdot a_i' + c\cdot a_j'
\end{aligned}$$

### 4.4 Gradient Computation

The gradient with respect to θ (using ∂c/∂θ = -s, ∂s/∂θ = c) is:

$$\frac{\partial\mathcal{L}}{\partial\theta} = 2\cdot\text{Re}\left[ (\partial a_i/\partial\theta)^* \cdot b_i + (\partial a_j/\partial\theta)^* \cdot b_j \right]$$

with

$$\begin{aligned}
\frac{\partial a_i}{\partial\theta} &= -s\cdot a_i' + e^{-i\phi}c\cdot a_j' \\
\frac{\partial a_j}{\partial\theta} &= -e^{i\phi}c\cdot a_i' - s\cdot a_j'
\end{aligned}$$

### 4.5 Adjoint Update

After computing gradients, we update adjoint entries for the next (earlier) layer by applying G†:

$$\begin{aligned}
b_i' &= c\cdot b_i + e^{-i\phi}s\cdot b_j \\
b_j' &= -e^{i\phi}s\cdot b_i + c\cdot b_j
\end{aligned}$$

All operations are O(1) per rotation. Total backward pass: O(M).

## 5. Empirical verification of the birth tree
The birth tree is not only a theoretical construct: it is observed in practice when tracking support evolution in sparse Givens simulations. In the measured examples below, the number of birth records equals $|\text{supp}(\psi_M)| - 1$ for final support, and the recorded structure matches the exact support.

| Mode | n | M | final_support | birth_count | birth_bytes | memory |
|------|---|----|---------------|-------------|-------------|--------|
| givens | 3 | 5 | 5 | 4 | 568 | 97 MB |
| decompose | 3 | 9 | 8 | 7 | 960 | 97 MB |

The table confirms that `birth_bytes` grows with the number of births. For RCS with larger $n$, this can become a bottleneck, but the trajectory method provides a constant-memory alternative for sampling.

## 6. Complexity Summary

| Phase | Cost |
|-------|------|
| Forward pass | O(M) |
| Expectation (L observables) | O(L·M) |
| Adjoint initialization | O(L·M) |
| Backward pass | O(M) |
| Memory | O(M) |

**No exponential anywhere.**

## 7. Limitations

The method does not apply to:
- Circuits containing non-Givens gates (CNOT, Hadamard, etc.)
- Circuits requiring O(2^n) Givens rotations to decompose
- General quantum algorithms like Shor's or random circuit sampling

In particular, our current `decompose` experiments show that a random Sycamore-like circuit with `n=10` and `depth=4` already reaches `final_support = 1024`, i.e. full support. This confirms that random two-qubit gates can destroy the sparsity advantage when they are decomposed into global Givens rotations.

However, for the specific class of Givens-only circuits, the simulation is both exact and efficient.

## 8. Conclusion

The principal result is not merely that the simulation runs in $O(M)$ time. It is that the simulator can record a **birth tree** of newly created amplitudes and thereby recover the structure of zeros and non-zeros in the state.

The birth tree makes support structure explicit: each new non-zero amplitude has a unique parent, and the directed forest of births can be reconstructed from the rotation parameters. This opens a path to compressed representations of quantum states beyond naive sparse vectors.

For circuits with small support growth, such as UCCSD and global Givens schemes, the tree is trivial and exact reconstruction is efficient. For RCS, the tree becomes exponential, but the trajectory method allows approximate sampling with constant memory.

This work therefore identifies a new structural lens on quantum simulation: exact support recovery via a causality graph, not just fast arithmetic.