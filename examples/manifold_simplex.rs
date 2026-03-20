//! Fisher-Rao simplex as a Riemannian manifold via `skel::Manifold`.
//!
//! Demonstrates that the probability simplex under the Fisher-Rao metric
//! satisfies the standard Riemannian manifold identities: exp/log round-trip,
//! geodesic distance from log norm, and parallel transport norm preservation.
//!
//! ```cargo
//! [dependencies]
//! infogeom = { path = "..", features = ["manifold"] }
//! skel = "0.1"
//! ndarray = "0.16"
//! ```

use infogeom::FisherRaoSimplex;
use ndarray::array;
use skel::Manifold;

fn main() {
    let m = FisherRaoSimplex::default();

    let p = array![0.7, 0.2, 0.1];
    let q = array![0.1, 0.2, 0.7];

    // log_map: find tangent vector from p to q
    let v = m.log_map(&p.view(), &q.view());
    println!("p = {p}");
    println!("q = {q}");
    println!("log_p(q) = {v:.6}");

    // exp_map: recover q from p + v
    let q_recovered = m.exp_map(&p.view(), &v.view());
    let err: f64 = q.iter()
        .zip(q_recovered.iter())
        .map(|(a, b)| (a - b).abs())
        .sum();
    println!("exp_p(log_p(q)) = {q_recovered:.6}");
    println!("round-trip L1 error: {err:.2e}");

    // Fisher-Rao distance via infogeom
    let d_rao = infogeom::rao_distance_categorical(
        p.as_slice().unwrap(),
        q.as_slice().unwrap(),
        1e-12,
    )
    .unwrap();
    println!("\nFisher-Rao distance: {d_rao:.6} rad");

    // Geodesic interpolation at t=0.25, 0.5, 0.75
    println!("\ngeodesic interpolation:");
    for &t in &[0.25, 0.5, 0.75] {
        let vt = v.mapv(|vi| vi * t);
        let pt = m.exp_map(&p.view(), &vt.view());
        println!("  t={t:.2}: {pt:.4}");
    }

    // Parallel transport: move v along the geodesic from p to q
    let tv = m.parallel_transport(&p.view(), &q.view(), &v.view());
    println!("\nparallel transport of v from T_p to T_q:");
    println!("  v  at p: {v:.6}");
    println!("  Pv at q: {tv:.6}");

    // Verify norm preservation (in the Fisher metric)
    let norm_at_p: f64 = p
        .iter()
        .zip(v.iter())
        .map(|(&pi, &vi)| vi * vi / (4.0 * pi))
        .sum::<f64>()
        .sqrt();
    let norm_at_q: f64 = q
        .iter()
        .zip(tv.iter())
        .map(|(&qi, &tvi)| tvi * tvi / (4.0 * qi))
        .sum::<f64>()
        .sqrt();
    println!("  Fisher norm at p: {norm_at_p:.6}");
    println!("  Fisher norm at q: {norm_at_q:.6}");
    println!("  norm preserved: {}", (norm_at_p - norm_at_q).abs() < 1e-6);
}
