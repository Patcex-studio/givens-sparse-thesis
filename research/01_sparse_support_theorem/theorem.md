# Theorem 1 — Sparse Support Bound

## Statement
For any sequence of `M` complex Givens rotations applied to `|0...0⟩`:

```
|supp(ψ_k)| ≤ 1 + k    for all k = 0,…,M
```

## Proof (induction)

**Base:** `k=0`, `ψ_0 = |0⟩` → `|supp| = 1`

**Step:** Assume true for `k`. Consider rotation `G(i,j,θ,φ)`.

- States `x ∉ {i,j}` unchanged → support preserved.
- If neither `i` nor `j` in support → neither becomes non‑zero (they remain zero); if exactly one is non‑zero, at most one new basis state appears.
- If one already present → support unchanged (+0).
- If both present → support unchanged (+0).

Thus:
```
|supp(ψ_{k+1})| ≤ |supp(ψ_k)| + 1 ≤ 1 + k + 1 = 1 + (k+1)
```

QED.

## Corollary
If `M = O(n)` then `|supp| = O(n)` everywhere.  
The state fits in a sparse map of polynomial size.