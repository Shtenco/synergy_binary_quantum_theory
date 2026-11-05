# BQG: held-out restricted-operator norm gate (2026-11-05)

## Exact scope

One real K5 Peter–Weyl state basis (source node 1, oriented triple (0,2,3), auxiliary channel 000, Jmax=7/2, one fixed logical column). Already known: 20 of 120 source basis inputs computed, G20 rank20, 1052 full outgoing keys, ||A P20 ψ||=0.06524151188803094, ||Q100 ψ||=0.08986860634523096.

At the current stage, no *upper* bound for ||A Q100|| exists. Claiming that 51% input weight is 51% output accuracy would be false.

## Rigorous one-column falsification theorem

Let the 120 source keys form an orthonormal finite-cutoff Hilbert basis `e_j`, and let Q project orthogonally onto the 100 uncomputed directions. If r is an index in Q, then Q e_r=e_r and ||e_r||=1; hence

```
||A Q|| = sup_{||v||=1} ||A Q v|| >= ||A Q e_r|| = ||A e_r||.
```

For actual source coefficient a_r !=0, the separately computed weighted full outgoing vector w_r=A(a_r e_r) has

```
||A e_r|| = ||w_r|| / |a_r|.
```

Therefore, if any independently held-out full real micro column obeys ||w_r||/|a_r| > 0.07259655461597475, no correct *operator norm* certificate of the form ||A Q||<=0.07259655 can exist. This does **not** imply ||A Q ψ|| is large, because the unknown 100 columns could cancel coherently on that one specific input.

The rigorous direct goal remains either actual ||A Q ψ|| from a complete 100-key output, or a dedicated upper bound on A restricted to the one-dimensional ray along Q ψ (rather than the entire Q100 space). An approximate spectral norm of the 20 already computed columns provides a lower bound on the global operator; it provides no upper bound on the unknown Q columns.

## Actual experiment

Workflow `bqg-node1-restricted-one-gram21.yml` selects the highest-weight held-out frozen KB coefficient with no post-hoc fitting, performs the original full third C(K) operator on that **one** basis state using `bqg_lorentzian_ka_partition.py`, verifies source SHA256, complete outgoing labels/spins, dividing by the actual coefficient, then calculates a genuine Gram G21 if the job completes. Logs include full NPZ in SHA-checked base64 to bypass artifact quota; source branch main is untouched.

No physical projector, independent HDA, interacting TT pole, or post-summation charged-basis closure claimed. A failure to compute or numerical cutoff failure must be recorded explicitly, not labelled a physics no-go.

## Expected decision tree

- If held-out ||Ae_r||>0.07259655, the sufficient criterion ||AQ||<=threshold is mathematically impossible. The next strategy must be **coherent remainder**, not a harder global spectrum search.
- If held-out ||Ae_r||<=threshold, the 10% method remains unverified, not proven; another e_r or a complete Q-restricted norm upper certificate is required.
- If direct task times out, profile that one costly SU2 source and persist exact intermediate recoupling data; do not alter Jmax or truncate high spins to get an apparently successful result.

Scientific status: full Lorentzian HDA, physical rigging projector, connected massless interacting TT pole are all OPEN.
