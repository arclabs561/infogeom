//! Geodesics and distances on the probability simplex.
//!
//! Demonstrates the three geodesic families (Fisher-Rao, mixture, exponential),
//! Rao and Hellinger distances, Fisher information diagonal, and natural gradient.

use infogeom::{
    e_geodesic, fisher_information_diagonal, fisher_rao_geodesic, hellinger, m_geodesic,
    natural_gradient, rao_distance_categorical,
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

    // --- Fisher information diagonal ---
    // For categorical distributions, the Fisher information matrix is diagonal:
    // F_ii = 1 / p_i. Categories with small probability have high Fisher information,
    // meaning the likelihood is more sensitive to changes in those coordinates.
    let fisher = fisher_information_diagonal(&p, tol).unwrap();
    println!("Fisher information diagonal (1/p_i):");
    for (i, (&pi, &fi)) in p.iter().zip(fisher.iter()).enumerate() {
        println!("  category {i}: p={pi:.2}  F_ii={fi:.4}");
    }
    println!("  Rare categories (small p_i) carry more Fisher information.");
    println!();

    // --- Natural gradient ---
    // The natural gradient rescales the Euclidean gradient by the inverse Fisher metric.
    // For categoricals, F^{-1} = diag(p_i), so natural_grad_i = p_i * grad_i.
    // This removes the implicit Euclidean geometry and follows the steepest direction
    // on the Fisher-Rao manifold. In practice, the natural gradient downweights
    // updates to rare categories and upweights updates to common ones.
    let grad = [1.0, -0.5, 0.3];
    let ng = natural_gradient(&p, &grad).unwrap();
    println!("Euclidean gradient: {grad:?}");
    println!("Natural gradient:   {ng:?}");
    println!("  natural_grad_i = p_i * grad_i (inverse Fisher metric applied)");
    println!(
        "  Category 0 (p=0.70): gradient 1.0 -> {:.4} (upweighted by large p)",
        ng[0]
    );
    println!(
        "  Category 2 (p=0.10): gradient 0.3 -> {:.4} (downweighted by small p)",
        ng[2]
    );

    println!();
    println!("All checks passed.");
}
