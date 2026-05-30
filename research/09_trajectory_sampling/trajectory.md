# Trajectory Sampling for Givens Circuits

## Definition

Trajectory sampling is an approximate simulation method for Givens circuits in which a single basis state is maintained at each step, and the state is collapsed after every rotation. The method is useful when only samples are required and not the full wavefunction.

## Pseudocode

```python
def sample_trajectory(circuit, n_qubits):
    state = {0: 1.0}
    for (gate_type, *params) in circuit:
        if gate_type == 'global_givens':
            i, j, theta, phi = params
            # apply the rotation to the current basis amplitude
            # compute probabilities for outcomes i and j
            # collapse to one of the basis states
        else:
            # apply the decomposed local Givens rotation
            # collapse within the affected 2-qubit block
    return bitstring
```

A more concrete implementation for a direct global Givens step is:

```python
def sample_trajectory_step(index, amplitude, i, j, theta, phi):
    # determine the two basis indices that participate in this rotation
    # calculate their amplitudes after the rotation
    # choose one of the two basis states randomly according to probability
    # return the new basis index with amplitude 1
```
```

## Complexity

- Time per sample: $O(M)$.
- Memory: $O(1)$.
- Total time for $N_{\text{samples}}$ samples: $O(M\cdot N_{\text{samples}})$.

## Applicability

The trajectory method is suitable for problems where samples are the target output rather than the full quantum state. It is particularly relevant for:

- random circuit sampling (RCS)
- quantum supremacy-style verification
- heuristic estimation of output distributions when exact state reconstruction is infeasible

Exact amplitude recovery still requires the full birth tree or the full sparse state, but trajectory sampling provides a constant-memory approximation for sampling.
