# Manifesto: The Givens Loophole

## They said:
> "Quantum states live in 2ⁿ dimensions. Classical simulation is impossible beyond 50 qubits."

## We noticed:
> "But what if the state was built only from rotations that touch two basis states at a time?"

## The blind spot:
Everyone assumed that because *some* quantum circuits fill the entire Hilbert space, *all* circuits must eventually do so.  
No one asked: *"How fast can the support grow?"*

The answer: linearly.  
`|supp(ψ_k)| ≤ 1 + k`

That's it. That's the whole trick.

## Why it's not trivial:
- It's exact (no approximations, no tensor networks, no Monte Carlo).
- It's deterministic.
- It runs in linear time.
- It handles backpropagation (gradients) in linear time.

## Why they missed it:
Because they were busy compressing the state. We just stopped storing zeros.

## This repository is not a code dump.  
## It is a mathematical proof that a non‑trivial class of quantum circuits is classically trivial.