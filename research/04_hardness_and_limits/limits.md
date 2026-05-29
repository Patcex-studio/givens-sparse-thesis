# Where the Method Breaks

## 1. Non‑Givens gates
If a single gate touches 4 basis states (e.g., CNOT, Hadamard on two qubits), the bound fails.  
One can construct circuits with `O(n)` such gates whose support becomes `2^{Ω(n)}`.

## 2. Decomposition overhead
Any unitary can be decomposed into `O(n²)` Givens rotations.  
If the original circuit has `O(n)` arbitrary gates, the Givens count becomes `O(n²)`, so support `O(n²)` — still polynomial, but not linear.

## 3. Intermediate explosion
The theorem bounds *final* support. It also bounds intermediate support (by the same induction).  
So no intermediate explosion **for Givens circuits**. This is a feature, not a limitation.

## 4. What remains exponential
- Shor's algorithm (requires exponentially many Givens rotations in any decomposition).
- Random circuit sampling (Google Sycamore) — each two‑qubit gate is not a Givens rotation; after decomposition, the number of Givens rotations is exponential in depth.

## Summary
The method trivializes **Givens‑only circuits** with `O(n)` rotations.  
It does **not** trivialize all of BQP.