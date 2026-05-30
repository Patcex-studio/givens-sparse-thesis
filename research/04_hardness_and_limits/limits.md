# Where the Method Breaks

## 1. Non-Givens gates
If a single gate touches 4 basis states (e.g., CNOT, Hadamard on two qubits), the bound fails.  
One can construct circuits with `O(n)` such gates whose support becomes `2^{\Omega(n)}`.

## 2. Decomposition overhead
Any unitary can be decomposed into `O(n^2)` Givens rotations.  
If the original circuit has `O(n)` arbitrary gates, the Givens count becomes `O(n^2)`, so support `O(n^2)` -- still polynomial, but not linear.

## 3. Global Givens versus two-qubit blocks
The support-bound proof applies to a sequence of Givens rotations taken as the original circuit. It does not apply to an arbitrary exact decomposition of a generic two-qubit gate into global Givens rotations.

### 3.1 Formal definitions
- **Computational basis.** For $n$ qubits we index basis states by binary strings $x \in \{0,1\}^n$. We identify a string with the integer $x \in [0,2^n)$.

- **Global Givens rotation.** A rotation $G(i,j,\theta,\phi)$ is the $2^n \times 2^n$ unitary that acts as the identity on every basis state except $|i\rangle$ and $|j\rangle$, with $i \neq j$. Its matrix is

$$
G(i,j,\theta,\phi) = I_{2^n} - (1-\cos\theta)(|i\rangle\langle i| + |j\rangle\langle j|) - \sin\theta\,(e^{-i\phi}|i\rangle\langle j| + e^{i\phi}|j\rangle\langle i|).
$$

Equivalently, the matrix has a single non-trivial $2 \times 2$ block at rows and columns $i,j$, and all other entries are $\delta_{ab}$.

- **Two-qubit gate on qubits $a,b$.** Let $U_{ab}$ be a unitary acting on qubits $a,b$ while leaving the remaining $n-2$ qubits untouched. In the full $2^n \times 2^n$ matrix it is block-diagonal: for each pattern $y \in \{0,1\}^{n-2}$ of the other qubits there is a $2 \times 2$ block $U_{ab}$. Hence there are exactly

$$
B = 2^{n-2}
$$

identical blocks.

- **Random circuit sampling (RCS).** An RCS circuit consists of depth-$d$ layers; each layer contains a random single-qubit SU(2) gate on every qubit and a random two-qubit SU(4) gate on each of the $n-1$ nearest-neighbour pairs. The two-qubit gates are Haar-random, so with probability 1 they are generic: the associated $2 \times 2$ block is not diagonalizable by a single Givens rotation and is not a scalar multiple of the identity.

### 3.2 Lemma 1 - A global Givens rotation touches exactly one block
Consider a global Givens rotation $G(i,j,\theta,\phi)$. Let

$$
S(i) = \{x \in \{0,1\}^n \mid x \text{ agrees with } i \text{ on all bits except possibly those where } i \text{ and } j \text{ differ}\}.
$$

The indices $i$ and $j$ differ on at least one bit. The bits where they agree define a pattern $p \in \{0,1\}^{n-2}$ for the untouched qubits. All basis states affected by the rotation share the same pattern $p$ on those other qubits. Any state with a different pattern lies in a different block of a two-qubit gate acting on $a,b$. Consequently the non-trivial $2 \times 2$ block of $G$ is confined to a single one of the $B$ blocks of that two-qubit gate.

*Proof.* The matrix of $G$ differs from the identity only in rows and columns $i$ and $j$. For any other basis index $x \neq i,j$, the corresponding column and row are unchanged, so $G$ acts trivially on the subspace defined by fixing all bits outside the positions where $i$ and $j$ differ. Therefore the support of $G$ is a single pair of basis states, i.e. a single block. QED

### 3.3 Lemma 2 - Lower bound on global Givens decomposition of a generic two-qubit gate
Let $U_{ab}$ be a generic unitary acting on qubits $a,b$ (a Haar-random element of SU(4)). Any exact decomposition

$$
U_{ab} = G_M G_{M-1} \cdots G_1,
$$

where each $G_k$ is a global Givens rotation, satisfies

$$
M \ge B = 2^{n-2}.
$$

*Proof.* Write the block-diagonal form of $U_{ab}$ as

$$
U_{ab} = \bigoplus_{y\in\{0,1\}^{n-2}} U_{ab}.
$$

Each direct summand corresponds to a fixed pattern $y$ of the remaining qubits, and every block is the same $2 \times 2$ matrix $U_{ab}$. By Lemma 1, each global Givens rotation affects at most one of the $B$ blocks. If $M < B$, at least one block would remain untouched, and the product would equal the identity on that block. This contradicts the generic non-triviality of $U_{ab}$ on every block. Hence $M \ge B$. 

*Remark.* The bound is essentially tight: one can implement $U_{ab}$ by applying the same two-qubit unitary to each block, i.e. by using $B$ global Givens rotations that differ only in the pattern of the untouched qubits.

## 4. Consequences for random circuit sampling
A depth-$d$ RCS circuit contains $d(n-1)$ nearest-neighbour two-qubit gates. By Lemma 2 each such gate requires at least $2^{n-2}$ global Givens rotations in any exact decomposition. Therefore the total number of rotations needed to represent the whole circuit satisfies

$$
M \ge d(n-1) 2^{n-2} = \Omega(2^n).
$$

This violates the constraint $M \le \text{poly}(n,d)$ for every nontrivial $n$ in the hard RCS regime. Hence no exact polynomial-size representation with global Givens rotations exists for generic RCS.

### 4.1 Empirical RCS observation
The current simulator confirms this behaviour. For a Sycamore-like decomposition with `n=10`, both `depth=4` and `depth=10` produce `final_support = 1024`, i.e. full support on $2^{10}$ amplitudes.

## 5. Why the earlier support-bound argument does not apply
The support-bound proof assumes a sequence of Givens rotations that are local to a fixed pair of qubits. In that setting the two states mixed by each rotation differ only on those qubits, and the other qubits remain fixed throughout the circuit. Then each rotation adds at most one new basis state to support.

In the present RCS decomposition problem, a global Givens rotation may mix two basis states that differ on any subset of qubits. A two-qubit gate is block-diagonal with respect to the remaining $n-2$ qubits, and by Lemma 1 each global Givens can affect only one block. The cost of reproducing every block is therefore exponential, so the support bound is irrelevant to the question of exact polynomial decomposition.

## 6. Summary
The method trivializes **Givens-only circuits** with `O(n)` rotations. It does **not** trivialize random circuit sampling or arbitrary BQP circuits.

For RCS, the exact decomposition into global Givens rotations requires at least $\Omega(2^n)$ rotations, because a generic two-qubit gate has $2^{n-2}$ identical blocks and each global Givens can touch only one of them.
