# Givens Sparse Thesis

**Claim:**  
Any quantum circuit built from complex Givens rotations admits exact classical simulation in `O(M)` time and memory, where `M` is the number of rotations.

**Core insight:**  
Each rotation touches exactly two basis states. Start with one non‑zero amplitude. After `M` rotations, you can never have more than `1 + M` non‑zero amplitudes. No explosion. No exponential. No approximation.

**This repository contains the formal development — no code, only ideas, proofs, and implications.**

---

## Navigation

- [Manifesto](./manifesto.md) — why this matters and why everyone missed it
- [Research notes](./research/01_sparse_support_theorem/theorem.md) — the support bound theorem
- [Algorithm](./research/02_linear_time_algorithm/algorithm.md) — linear‑time forward/backward simulation
- [Limits](./research/04_hardness_and_limits/limits.md) — where the method breaks
- [Synthesis](./synthesis/article.md) — formatted as a short paper

---

**License:** BSD 3-Clause  
**Author:** Jocer Speis 
**Status:** Theoretical development complete. Experimental validation in separate repository.