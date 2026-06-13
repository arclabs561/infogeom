# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
