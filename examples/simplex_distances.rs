//! Geodesics and distances on the probability simplex.
//!
//! Demonstrates the three geodesic families (Fisher-Rao, mixture, exponential)
//! alongside Rao and Hellinger distances.

use infogeom::{
    e_geodesic, fisher_rao_geodesic, hellinger, m_geodesic, natural_gradient,
    rao_distance_categorical,
};

fn main() {
    let p = [0.70, 0.20, 0.10];
    let q = [0.10, 0.20, 0.70];
    let tol = 1e-12;

    // --- Distances ---
    let d_rao = rao_distance_categorical(&p, &q, tol).unwrap();
    let d_hel = hellinger(&p, &q, tol).unwrap();

    println!("Rao distance (radians): {d_rao:.6}");
    println!("Hellinger distance:     {d_hel:.6}");
    println!();

    // --- Geodesic paths ---
    println!(
        "{:>4}  {:>28}  {:>28}  {:>28}",
        "t", "Fisher-Rao", "Mixture (a=-1)", "Exponential (a=+1)"
    );
    println!("{}", "-".repeat(96));

    for step in 0..=10 {
        let t = step as f64 / 10.0;
        let fr = fisher_rao_geodesic(&p, &q, t, tol).unwrap();
        let m = m_geodesic(&p, &q, t, tol).unwrap();
        let e = e_geodesic(&p, &q, t, tol).unwrap();

        println!(
            "{t:>4.1}  [{:.4}, {:.4}, {:.4}]  [{:.4}, {:.4}, {:.4}]  [{:.4}, {:.4}, {:.4}]",
            fr[0], fr[1], fr[2], m[0], m[1], m[2], e[0], e[1], e[2],
        );
    }

    println!();

    // --- Geodesic midpoint verification ---
    let mid = fisher_rao_geodesic(&p, &q, 0.5, tol).unwrap();
    let d_pm = rao_distance_categorical(&p, &mid, tol).unwrap();
    let d_mq = rao_distance_categorical(&mid, &q, tol).unwrap();
    println!("Fisher-Rao midpoint check:");
    println!("  d(p, mid) = {d_pm:.6}");
    println!("  d(mid, q) = {d_mq:.6}");
    println!("  d(p, q)   = {d_rao:.6}");
    assert!((d_pm + d_mq - d_rao).abs() < 1e-10);
    println!();

    // --- Natural gradient ---
    let grad = [1.0, -0.5, 0.3];
    let ng = natural_gradient(&p, &grad).unwrap();
    println!("Euclidean gradient: {grad:?}");
    println!("Natural gradient:   {ng:?}");
    println!("  (each component scaled by p_i)");

    println!();
    println!("All checks passed.");
}
