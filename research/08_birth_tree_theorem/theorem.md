# Birth Tree Theorem

**Theorem (Birth Tree Theorem).**

For any sequence of $M$ complex Givens rotations applied to the initial state $|0\ldots0\rangle$, the process of appearance of new non-zero amplitudes defines a directed tree (forest) such that:

- Roots are basis states that never appear as a child of another birth event (typically a single root $|0\ldots0\rangle$).
- A directed edge $(\text{parent} \to \text{child})$ corresponds to a rotation at step $t$ in which the parent was already non-zero and the child became non-zero.
- The number of vertices in the tree is $|\text{supp}(\psi_M)|$.
- The number of edges is $|\text{supp}(\psi_M)| - 1$.
- Each edge can be recovered from the rotation parameters and the amplitude of the child (inverse transformation).

## Proof

The proof proceeds by establishing that each new non-zero amplitude can only arise from a single existing non-zero amplitude.

1. Assume a new amplitude appears at step $t$ as a result of applying a Givens rotation $G(i,j,\theta,\phi)$. The rotation acts on exactly two basis indices $i$ and $j$.
2. At most one of the two amplitudes $a_i$ and $a_j$ can be zero before the rotation while the other is non-zero. If both were zero, the rotation cannot create a new non-zero amplitude; if both were non-zero, the amplitude is not new.
3. Therefore, the new non-zero amplitude at step $t$ is created from exactly one existing non-zero amplitude. That existing amplitude is the unique parent of the new child.
4. Since step indices increase monotonically, no cycle can form: a child can only be born at a later step than its parent.
5. Consequently, the birth relations form a directed acyclic graph in which every vertex except roots has exactly one parent, i.e., a directed forest. Because the initial state has one non-zero amplitude and each birth introduces one new vertex, the graph has $|\text{supp}(\psi_M)|$ vertices and $|\text{supp}(\psi_M)| - 1$ edges.

Finally, each birth record stores the rotation parameters $(i,j,\theta,\phi)$ and the child amplitude after the step. The inverse Givens transformation reconstructs the parent amplitude from the child amplitude when $\sin\theta \neq 0$, so every edge is recoverable.

## Consequences

- The memory required to store the birth tree is $O(|\text{supp}|)$ in the worst case.
- The full amplitude structure at all steps can be recovered in $O(M)$ operations given only the birth tree and the rotation parameters.
- For circuits with small support growth (e.g., UCCSD, global Givens schemes), the birth tree is trivial.
- For circuits with $|\text{supp}| = 2^n$ (e.g., RCS), the birth tree is exponential, but it may be bypassed by the trajectory method described in `research/09_trajectory_sampling/trajectory.md`.

*Note.* The birth tree argument does not imply that every two-qubit gate can be decomposed into a polynomial number of global Givens rotations. A generic two-qubit gate has $2^{n-2}$ identical blocks, and each global Givens rotation can affect only one block, so exact decomposition generally requires exponential resources.

## Empirical validation

Measured experiments using `venv\Scripts\python.exe rcs_givens.py --track-births` confirm the birth-tree property for small random circuits: each final non-zero amplitude has a unique birth parent, and the number of recorded births equals $|\text{supp}(\psi_M)| - 1$ in the tracked examples.
