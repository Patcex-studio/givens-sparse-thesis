# FAQ

## How is this different from tensor network methods?

Tensor networks (e.g., MPS, PEPS) compress the quantum state by exploiting low entanglement. They provide approximate simulations with error bounds. Our method is **exact** — no approximation, no truncation. The sparsity comes from the circuit structure (Givens rotations), not from entanglement assumptions.

## Why doesn't this refute quantum supremacy?

Quantum supremacy experiments (e.g., Google Sycamore) use random circuits with arbitrary two-qubit gates. These gates are **not** Givens rotations. Decomposing a generic two-qubit gate into Givens rotations requires `O(n²)` rotations, and the support bound becomes `O(n²)`, not linear. Moreover, the decomposition itself is exponential in circuit depth for random circuits. Our result applies only to circuits composed **solely** of Givens rotations.

## Why does `decompose` on `n=10`, `depth=4` still reach full support?

This is exactly the point. In the current implementation, a Sycamore-like decomposition on `n=10` with `depth=4` produces `final_support = 1024`. That means the circuit has already explored the full $2^{10}$ basis within the decomposed Givens representation, so the sparse support argument does not provide any compression in this regime.

The reason is that the original support bound holds for a sequence of local Givens rotations taken as the actual circuit. A generic two-qubit gate decomposed into global Givens rotations is a different problem: each global rotation can affect only one block among $2^{n-2}$ blocks, and a generic block-diagonal gate requires exponentially many such rotations.

## Can I use this to simulate Shor's algorithm?

No. Shor's algorithm uses a variety of gates including Hadamards, CNOTs, and controlled-phase gates. While each can be decomposed into Givens rotations, the number of rotations required is exponential in the number of qubits. The support bound then becomes exponential, and the linear-time advantage is lost.

## You mention "experimental validation in a separate repository" — where is the code?

The theoretical development in this repository is code-agnostic. The implementation (C++/CUDA/Python) is in a separate repository: `givens-sparse-impl`. This keeps the research notes clean and focused on mathematics.

## Is the method applicable to linear optics quantum computing (LOQC)?

Yes, with caveats. LOQC uses beam splitters and phase shifters, which are essentially Givens rotations. However, photon loss and detection inefficiency introduce noise that our exact simulation does not account for. The method applies to the ideal, lossless LOQC model.

## What about circuits with a mix of Givens and non-Givens gates?

If a circuit contains `k` non-Givens gates, each can be decomposed into `O(n²)` Givens rotations. The total number of Givens rotations becomes `O(n²·depth)`, and the support bound is `O(n²·depth)`. This is still polynomial, but not linear. The linear-time advantage holds only when the original circuit has `O(n)` Givens rotations.

## Why complex Givens rotations and not just real ones?

Real Givens rotations (θ only) can only generate a subset of unitaries. Complex Givens rotations (θ, φ) provide full control over the two-dimensional subspace, enabling universal quantum computation when combined with single-qubit gates. The theorem holds for the more general complex case.

## Does this work for other parameterized quantum circuits?

The key property is that each gate affects a small, fixed number of basis states. Any gate that touches `d` basis states would lead to a support bound of `O(d·M)`. Givens rotations are special because `d = 2`. Gates affecting more basis states (e.g., CNOT affects 4 states) do not yield the same linear bound.

## What's next?

Potential directions:
- Extend to other structured gate sets (e.g., Clifford gates)
- Hybrid algorithms combining sparse simulation with tensor networks
- Automatic detection of Givens-only subcircuits in larger algorithms