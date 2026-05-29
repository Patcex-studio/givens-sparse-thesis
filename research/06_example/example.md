# Numerical Example: 3 Givens Rotations on 2 Qubits

This example demonstrates the sparse support theorem and linear-time algorithm with a concrete numerical example.

## Setup

- **Number of qubits**: `n = 2`
- **Initial state**: `|ψ₀⟩ = |00⟩`
- **Number of rotations**: `M = 3`

## Rotation Sequence

We apply three Givens rotations:

| Step | i | j | θ | φ |
|------|---|---|---|---|
| 1 | 0 | 1 | π/4 | 0 |
| 2 | 0 | 1 | π/4 | π/2 |
| 3 | 0 | 1 | π/4 | 0 |

## Step-by-Step Execution

### Initial State
```
|ψ₀⟩ = |00⟩
Amplitudes: a_00 = 1, all others = 0
Support: {00} (size = 1)
```

### Rotation 1: G(0,1,π/4,0)

Parameters: `c = cos(π/4) = √2/2`, `s = sin(π/4) = √2/2`, `e = e^{-i·0} = 1`

Forward update:
```
a_00' = c·a_00 - s·e·a_01 = (√2/2)·1 - (√2/2)·0 = √2/2
a_01' = s·conj(e)·a_00 + c·a_01 = (√2/2)·1 + (√2/2)·0 = √2/2
```

State after rotation 1:
```
|ψ₁⟩ = (√2/2)|00⟩ + (√2/2)|01⟩
```

Support: `{00, 01}` (size = 2)

> Note: At most one new basis state is created per rotation, as guaranteed by the sparse support theorem.

### Rotation 2: G(0,1,π/4,π/2)

Parameters: `c = √2/2`, `s = √2/2`, `e = e^{-i·π/2} = -i`

Forward update:
```
a_00' = c·a_00 - s·e·a_01 = (√2/2)·(√2/2) - (√2/2)·(-i)·(√2/2)
      = 1/2 + i/2 = (1+i)/2

a_01' = s·conj(e)·a_00 + c·a_01 = (√2/2)·(i)·(√2/2) + (√2/2)·(√2/2)
      = i/2 + 1/2 = (1+i)/2
```

State after rotation 2:
```
|ψ₂⟩ = ((1+i)/2)|00⟩ + ((1+i)/2)|01⟩
```

Support: `{00, 01}` (size = 2, no new states added)

### Rotation 3: G(0,1,π/4,0)

Same as rotation 1: `c = √2/2`, `s = √2/2`, `e = 1`

Forward update:
```
a_00' = c·a_00 - s·e·a_01 = (√2/2)·((1+i)/2) - (√2/2)·1·((1+i)/2)
      = (√2/2)·(1+i)/2 - (√2/2)·(1+i)/2 = 0

a_01' = s·conj(e)·a_00 + c·a_01 = (√2/2)·1·((1+i)/2) + (√2/2)·((1+i)/2)
      = (√2/2)·(1+i) = (√2/2 + i√2/2)
```

State after rotation 3:
```
|ψ₃⟩ = 0·|00⟩ + ((√2/2)·(1+i))|01⟩
```

Support: `{01}` (size = 1)

## Final State

```
|ψ₃⟩ = 0·|00⟩ + (√2/2 + i√2/2)·|01⟩
```

The support size never exceeded `1 + M = 1 + 3 = 4`, and in fact stayed ≤ 2 throughout.

## Expectation Value Calculation

Compute `⟨Z₁⟩` (Pauli Z on qubit 1):

```
⟨Z₁⟩ = ⟨ψ₃| (I⊗Z) |ψ₃⟩
      = |a_00|²·⟨00|I⊗Z|00⟩ + |a_01|²·⟨01|I⊗Z|01⟩
        + cross terms (which vanish for Z)
```

Since `Z|0⟩ = |0⟩` and `Z|1⟩ = -|1⟩`:

```
⟨Z₁⟩ = |a_00|²·1 + |a_01|²·(-1)
      = 0² - |√2/2 + i√2/2|²
      = - |(√2/2)(1+i)|²
      = - (1/2)·|1+i|²
      = - (1/2)·2
      = -1
```

This matches the expectation: the final state is purely |01⟩, so qubit 1 is in state |1⟩ with Z-eigenvalue -1.

## Verification of Sparse Support Bound

At each step, the support size never exceeds `1 + k` (where `k` is the number of rotations applied so far):

- Step 0: support size = 1 ≤ 1 + 0 = 1 ✓
- Step 1: support size = 2 ≤ 1 + 1 = 2 ✓
- Step 2: support size = 2 ≤ 1 + 2 = 3 ✓
- Step 3: support size = 1 ≤ 1 + 3 = 4 ✓

The tight bound `|supp(ψ_k)| ≤ 1 + k` holds at every step.