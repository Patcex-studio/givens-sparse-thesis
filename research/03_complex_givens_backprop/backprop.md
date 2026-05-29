# Complex Givens Backpropagation

## Rotation matrix
```
G(θ,φ) = [ c      -e^{-iφ} s ]
         [ e^{iφ} s   c       ]
```
where `c = cosθ`, `s = sinθ`.

## Inverse (adjoint)
```
G†(θ,φ) = [ c      e^{-iφ} s ]
          [ -e^{iφ} s   c     ]
```

## Forward update
```
a_i' = c·a_i - e^{-iφ} s·a_j
a_j' = e^{iφ} s·a_i + c·a_j
```

## Backward (reverse pass)
Given post‑rotation `a_i', a_j'` and adjoint entries `b_i, b_j`:

Pre‑rotation amplitudes:
```
a_i =  c·a_i' + e^{-iφ} s·a_j'
a_j = -e^{iφ} s·a_i' + c·a_j'
```

Gradient w.r.t. `θ` (using `∂c/∂θ = -s`, `∂s/∂θ = c`):
```
∂L/∂θ = 2·Re[ (∂a_i/∂θ)*·b_i + (∂a_j/∂θ)*·b_j ]
```
with
```
∂a_i/∂θ = -s·a_i' + e^{-iφ} c·a_j'
∂a_j/∂θ = -e^{iφ} c·a_i' - s·a_j'
```

Update adjoint for next layer (apply `G†` to `b`):
```
b_i' =  c·b_i + e^{-iφ} s·b_j
b_j' = -e^{iφ} s·b_i + c·b_j
```

All operations `O(1)` per rotation.

> Note: The gradient can be expressed directly in terms of `a_i', a_j'` without explicitly computing `a_i, a_j`. The formulas above are provided for clarity; an implementation may compute `∂L/∂θ` using only `a_i', a_j', b_i, b_j` and the trigonometric functions.