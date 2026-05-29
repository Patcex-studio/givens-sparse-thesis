# Givens Sparse Thesis

[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-coming_soon-red.svg)](https://arxiv.org)
[![GitHub release](https://img.shields.io/github/v/release/Patcex-studio/givens-sparse-thesis)](https://github.com/Patcex-studio/givens-sparse-thesis/releases)
[![DOI](https://img.shields.io/badge/DOI-coming_soon-blue)](https://zenodo.org)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

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
## Commercial use

This software is released under the BSD 3-Clause license, which **permits commercial use**.

Pharmaceutical, biotech, and material science companies are welcome to:
- Integrate this method into their internal simulation pipelines
- Build proprietary products based on this method
- Contact patcex@proton.me for enterprise support, training, or custom development

---

## Quick start

```bash
git clone https://github.com/Patcex-studio/givens-sparse-thesis.git
cd givens-sparse-thesis
python verify_givens.py   # runs the verification script
```

---

<img width="1672" height="941" alt="l;rk3-2ec" src="https://github.com/user-attachments/assets/fcf74496-be29-4919-8cb3-98918dc3fb03" />

---


**License:** BSD 3-Clause  
**Author:** Jocer Speis 
**Status:** Theoretical development complete. Experimental validation in separate repository.


