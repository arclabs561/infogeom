# Examples

## Which example should I run?

| I want to... | Example |
|---|---|
| Compare common distances on the simplex | `simplex_distances` |
| See how alpha-geodesics change an interpolation | `alpha_geodesic_family` |
| Diagnose retrieval distribution shift | `retrieval_distribution_shift` |
| Compare Fisher-Rao, Hellinger, KL, JS, and TV | `divergence_geometry` |
| Compare information-geometric interpolation with OT interpolation | `ot_interpolation_geometry` |
| Use the simplex through `skel::Manifold` | `manifold_simplex` |

## Example descriptions

- `simplex_distances`: prints Fisher-Rao and Hellinger distances, geodesic paths, Fisher information, and natural-gradient scaling for a three-category distribution.
- `alpha_geodesic_family`: renders midpoint distributions for mixture, Fisher-Rao, and exponential geodesics. Useful when you want to see why the alpha parameter changes the shape of the path.
- `retrieval_distribution_shift`: treats ranker outputs as categorical distributions over result families. Shows Fisher-Rao drift, midpoint behavior, and a natural-gradient feedback step.
- `divergence_geometry`: compares `infogeom` distances with `logp` divergences and verifies the shared Bhattacharyya-coefficient relationship.
- `ot_interpolation_geometry`: compares Fisher-Rao and mixture paths with a transport-plan interpolation from `wass`.
- `manifold_simplex`: exposes the Fisher-Rao simplex as a `skel::Manifold`. Requires the `manifold` feature.

## Running

```sh
cargo run --example retrieval_distribution_shift
cargo run --example alpha_geodesic_family
cargo run --features manifold --example manifold_simplex
```
