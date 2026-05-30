# Linear‑Time Exact Simulator

## Data structure
- Hash map / dictionary: `basis index → complex amplitude`
- Only non‑zero amplitudes stored.

## Forward pass
```
initialize map[0] = 1.0
for each rotation (i, j, θ, φ):
    a_i = map.get(i, 0)
    a_j = map.get(j, 0)
    c = cosθ, s = sinθ, e = exp(-iφ)
    a_i' = c*a_i - s*e*a_j
    a_j' = s*conj(e)*a_i + c*a_j
    map[i] = a_i'
    map[j] = a_j'
```
Each step: `O(1)` dictionary operations. Total `O(M)`.

## 3.3 Birth tracking algorithm
The simulator can also record the **birth tree** of newly created non-zero amplitudes. A birth record is a tuple `(step, parent, child, i, j, theta, phi)`.

```python
births = []  # list of tuples (step, parent, child, i, j, theta, phi)

def apply_rotation_pair_with_tracking(state, i, j, theta, phi, step, births):
    old_i = state.get(i, 0j)
    old_j = state.get(j, 0j)
    c = cos(theta)
    s = sin(theta)
    e = exp(-1j * phi)

    new_i = c * old_i - s * e * old_j
    new_j = s * e.conjugate() * old_i + c * old_j

    if new_i != 0 and i not in state:
        births.append((step, j, i, i, j, theta, phi))
    if new_j != 0 and j not in state:
        births.append((step, i, j, i, j, theta, phi))

    if abs(new_i) < eps:
        state.pop(i, None)
    else:
        state[i] = new_i
    if abs(new_j) < eps:
        state.pop(j, None)
    else:
        state[j] = new_j
```

Each birth record stores the parent index, child index, gate indices `i,j`, and the rotation parameters `θ, φ`. The birth tree is built incrementally during the forward pass.

## 3.4 Reconstructing amplitudes from the birth tree
A birth record determines the unique parent that created a new amplitude. If the child amplitude after the rotation is known, the parent amplitude before the rotation can be reconstructed by inverting the Givens rotation.

For a birth record with sine `s = \\sin θ` and phase `φ`, the inverse relation is:

$$
a_{\text{parent}}^{\text{before}} = \frac{1}{s} e^{-i φ} a_{\text{child}}^{\text{after}} \quad (s \neq 0).
$$

If $s = 0$, the rotation does not create a birth and no new child is introduced. Since each birth has a unique parent, the recorded edges form a directed forest, and the full structure of non-zero amplitudes can be recovered in $O(M)$ time when combined with the rotation parameters.

> Note: this support bound applies to a sequence of Givens rotations taken as the original circuit. It does not imply that an arbitrary generic two-qubit gate can be decomposed into a polynomial number of global Givens rotations without losing the sparse structure.

## Expectation values
For Pauli string `P` with mask `m` and phase `z`:
```
⟨ψ|P|ψ⟩ = Σ_{x∈supp} ψ_x^* ψ_{x⊕m} (-1)^{x·z}·i^{|m∧z|}
```
Sum over `O(|supp|) = O(M)` terms.  
For `L` observables: `O(L·M)`.

## Backpropagation (gradients)

**Initialization of adjoint vector**  
For the last layer (`k = M`), the adjoint entries `b_x` are computed as:
```
b_x = 2 · (H ψ_M)_x
```
where `H` is the observable (a sum of Pauli strings). Since `ψ_M` is sparse (support size ≤ 1+M), and each Pauli string maps support to at most `|supp|` terms, the total cost is `O(L·M)`.

**Reverse pass invariant**  
At the start of processing rotation `k` (going backwards), `b_i, b_j` are the adjoint variables **after** the rotation has been applied to the state. The algorithm:
1. Computes the gradient `∂L/∂θ_k` using `b_i, b_j` and the post‑rotation amplitudes `a_i', a_j'` (stored from forward pass).
2. Updates `b` to obtain the adjoint variables **before** the rotation using `G†`.
These become the `b` for the next (earlier) rotation.

*Detailed derivation of gradient formulas is provided in `research/03_complex_givens_backprop/backprop.md`.*

Process rotations in reverse order. For each rotation `(i, j, θ, φ)` with post‑rotation amplitudes `a_i', a_j'` and adjoint entries `b_i, b_j`:

1. **Restore pre‑rotation amplitudes** (inverse Givens rotation):
   ```
   a_i =  c·a_i' + e^{-iφ} s·a_j'
   a_j = -e^{iφ} s·a_i' + c·a_j'
   ```

2. **Compute gradient** `∂L/∂θ` via chain rule (detailed derivation in `research/03_complex_givens_backprop/backprop.md`).

3. **Update adjoint entries** for the next (earlier) layer:
   ```
   b_i' =  c·b_i + e^{-iφ} s·b_j
   b_j' = -e^{iφ} s·b_i + c·b_j
   ```
   and store them back into the map.

Each step is `O(1)`. Total reverse pass: `O(M)`.

## Complexity summary
| Phase | Cost |
|-------|------|
| Forward | O(M) |
| Expectation (L observables) | O(L·M) |
| Adjoint init | O(L·M) |
| Backward | O(M) |
| Memory | O(M) |

**No exponential anywhere.**