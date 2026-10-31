# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Source history (unversioned)

These entries identify source changes without assigning a package release.

- Reject non-finite `alpha` in `alpha_geodesic` through the existing
  `Error::Logp(logp::Error::Domain("alpha must be finite"))` error, after
  validating the other inputs.
  ([a20e7f1](https://github.com/arclabs561/infogeom/commit/a20e7f142bd00049e9dc394c51907e4b5cf54d3f))

## [Unreleased]

### Fixed

- `alpha_geodesic` with `alpha = 1.0` and `tol = 0.0` now delegates to
  `e_geodesic`. Before, the strict comparison skipped the delegation, the
  embedding exponent became 0, and every point on the curve was uniform.

## [0.2.3] - 2026-10-30

### Fixed

- `fisher_orthogonal_project` and `fisher_orthogonal_project_full` scale the
  Cholesky ridge to the largest diagonal of `GᵀFG` instead of adding a fixed
  1e-12. With small gradients or a small diagonal Fisher (an EWC-style
  diagonal near 1e-8), the fixed ridge dominated and 50 to 99% of the past
  component survived projection.
- Categorical Fisher–Rao distance no longer snaps nearby distributions to zero
  using the validation tolerance, which could violate the triangle inequality.
  Accepted positive-mass inputs are normalized and evaluated with a square-root
  chord formula. Tolerance must be finite and nonnegative; zero total mass is
  rejected through the existing wrapped domain error. Geodesic, manifold and
  finite-alpha behavior are unchanged.

## [0.2.2] - 2026-10-09

### Changed
- Declare `rust-version = "1.75"`.

### Fixed
- `natural_gradient` now returns `p_i * (g_i - p.g)`, the Fisher-Rao gradient
  on the simplex tangent space. The previous `p_i * g_i` summed to `p.g` and
  stepped off the simplex. The result now sums to zero and ignores a constant
  shift of `g`, so outputs differ from 0.2.1 for every `g` that is not already
  centered under `p`.

## [0.2.1] - 2026-06-13

### Added
- `alpha_geodesic(p, q, t, alpha, tol)`: the general alpha-family geodesic on the simplex via the affine alpha-embedding, reducing exactly to `m_geodesic` at `alpha = -1` and `e_geodesic` at `alpha = +1`. `alpha > 1` requires strictly positive entries.
- `alpha_geodesic_family` example: renders the t=0.5 midpoint of two opposing humps under the mixture, Fisher-Rao, and exponential geometries.

## [0.2.0] - 2026-04-16

### Added
- Initial published release: information-geometry primitives.
- Fisher-orthogonal gradient projection.
- `skel::Manifold` implementation for the Fisher-Rao simplex.
- OT interpolation geometry example (with `wass`).
- `divergence_geometry` example.
- Property test for the Rao/Bhattacharyya formula.
- Quickstart docs and runnable examples.
