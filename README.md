# infogeom

[![crates.io](https://img.shields.io/crates/v/infogeom.svg)](https://crates.io/crates/infogeom)
[![Documentation](https://docs.rs/infogeom/badge.svg)](https://docs.rs/infogeom)
[![CI](https://github.com/arclabs561/infogeom/actions/workflows/ci.yml/badge.svg)](https://github.com/arclabs561/infogeom/actions/workflows/ci.yml)

Information geometry on the probability simplex.

## Quickstart

```toml
[dependencies]
infogeom = "0.2.1"
```

```rust
use infogeom::{fisher_rao_geodesic, rao_distance_categorical};

let p = [0.70, 0.20, 0.10];
let q = [0.10, 0.20, 0.70];

// Fisher-Rao geodesic midpoint.
let mid = fisher_rao_geodesic(&p, &q, 0.5, 1e-12).unwrap();

// The midpoint is equidistant from both endpoints.
let d_full = rao_distance_categorical(&p, &q, 1e-12).unwrap();
let d_half = rao_distance_categorical(&p, &mid, 1e-12).unwrap();
assert!((d_half - d_full / 2.0).abs() < 1e-10);
```

## API

### Distances

| Function | Description |
|---|---|
| `rao_distance_categorical(p, q, tol)` | Fisher-Rao distance on the simplex (radians, range [0, pi]) |
| `hellinger(p, q, tol)` | Hellinger distance (re-exported from `logp`, range [0, 1]) |

### Geodesics

All take `(p, q, t, tol)` where `t` in [0, 1] interpolates from `p` to `q` (`alpha_geodesic` takes an extra `alpha` argument).

| Function | Alpha | Description |
|---|---|---|
| `fisher_rao_geodesic` | 0 | Riemannian geodesic via sphere embedding (slerp) |
| `m_geodesic` | -1 | Mixture geodesic: linear interpolation in probability space |
| `e_geodesic` | +1 | Exponential geodesic: linear interpolation in log space (requires strictly positive entries) |
| `alpha_geodesic` | any | The alpha-family geodesic for arbitrary `alpha`; reduces exactly to `m_geodesic` at -1 and `e_geodesic` at +1. `alpha > 1` requires strictly positive entries |

### Fisher information and natural gradient

| Function | Description |
|---|---|
| `fisher_information_diagonal(p, tol)` | Diagonal of the Fisher information matrix: `[1/p_1, ..., 1/p_n]` |
| `natural_gradient(p, euclidean_grad)` | Natural gradient: `p_i * g_i` (inverse Fisher metric applied to Euclidean gradient) |

## Tolerances

- Inputs are validated as simplex distributions (nonnegative, sum approximately 1) using `tol`.
- The `tol` parameter also controls degenerate-case snapping (e.g., BC near 1.0).

## References

- Amari & Nagaoka (2000), *Methods of Information Geometry*. Ch. 2-3 (Fisher metric, alpha-connections).
- Amari (1998), "Natural Gradient Works Efficiently in Learning", *Neural Computation* 10(2).
- Frank Nielsen's information geometry portal: [franknielsen.github.io/IG](https://franknielsen.github.io/IG/index.html)

## Examples

- `cargo run --example simplex_distances`: geodesics, distances, and natural gradient
- `cargo run --example divergence_geometry`: cross-crate comparison with `logp` divergences
- `cargo run --example retrieval_distribution_shift`: Fisher-Rao drift and natural-gradient feedback for ranker output distributions
- `cargo run --example ot_interpolation_geometry`: geodesics vs. optimal transport interpolation
- `cargo run --example manifold_simplex`: Fisher-Rao simplex as a `skel::Manifold` (requires `--features manifold`)
- `cargo run --example alpha_geodesic_family`: the alpha-family midpoint of two humps, rendered: mixture keeps both peaks, exponential pulls to the overlap, Fisher-Rao between

## License

MIT OR Apache-2.0
