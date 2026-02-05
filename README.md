# infogeom

Information geometry primitives.

`logp` provides divergence/entropy functionals on the simplex; `infogeom` builds **geometric**
structure on top (metrics, distances, and later projections/dual connections) without mixing in
application policy.

This crate is **classical** information geometry (simplex / exponential-family direction), not
quantum IG (see the separate `qig` crate in the broader workspace for density-matrix geometry).

## Quickstart

Add to your `Cargo.toml`:

```toml
[dependencies]
infogeom = "0.1"
```

Compute Rao (Fisher–Rao) and Hellinger distances on the simplex:

```rust
use infogeom::{rao_distance_categorical, hellinger};

let p = [0.70, 0.20, 0.10];
let q = [0.10, 0.20, 0.70];

let d_rao = rao_distance_categorical(&p, &q, 1e-12).unwrap();
let d_hel = hellinger(&p, &q, 1e-12).unwrap();

assert!(d_rao >= 0.0);
assert!((0.0..=1.0).contains(&d_hel));
```

## API tour

- `rao_distance_categorical(p, q, tol) -> Result<f64>`
  - Fisher–Rao (Rao) distance on the simplex via the sphere embedding:
    \(d_{FR}(p,q) = 2\arccos(\sum_i \sqrt{p_i q_i})\)
- `hellinger(p, q, tol) -> Result<f64>`
  - Hellinger distance (bounded metric on the simplex) via `logp`:
    \(H^2(p,q) = 1 - \sum_i \sqrt{p_i q_i}\)

## Invariants and tolerances

- Inputs are validated as simplex distributions (nonnegative, sum≈1) using `tol`.
- Distances are returned as:
  - Rao: radians in \([0, \pi]\)
  - Hellinger: \([0, 1]\)

Background reading:
- Frank Nielsen’s “Information geometry and divergences” portal: https://franknielsen.github.io/IG/index.html

## Examples

- `cargo run --example simplex_distances`

## Roadmap (near-term)

- Simplex geodesic/interpolation helpers (sphere embedding).
- \(\alpha\)-geometry scaffolding and dual-flat primitives (staying policy-free).

## License

MIT OR Apache-2.0

