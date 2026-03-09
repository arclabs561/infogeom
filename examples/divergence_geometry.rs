//! Cross-crate example: `logp` divergences alongside `infogeom` geometric distances.
//!
//! Shows how f-divergences from `logp` relate to the Fisher--Rao geodesic distance
//! and Hellinger metric provided by `infogeom`. All computations use discrete
//! distributions over 4 categories.

use infogeom::{hellinger, rao_distance_categorical};

fn main() {
    let tol = 1e-12;

    // Three distributions on a 4-simplex, ordered by "distance" from p.
    let p = [0.40, 0.30, 0.20, 0.10];
    let q_near = [0.35, 0.30, 0.25, 0.10]; // small perturbation
    let q_far = [0.10, 0.10, 0.30, 0.50]; // large shift

    let pairs: &[(&str, &[f64; 4], &[f64; 4])] =
        &[("p vs q_near", &p, &q_near), ("p vs q_far", &p, &q_far)];

    println!(
        "{:<14} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "pair", "Rao", "Hellinger", "KL(p||q)", "JS", "TV", "Bhatt_d"
    );
    println!("{}", "-".repeat(74));

    for &(label, a, b) in pairs {
        let rao = rao_distance_categorical(a, b, tol).unwrap();
        let hel = hellinger(a, b, tol).unwrap();

        // logp divergences computed directly.
        let kl = logp::kl_divergence(a, b, tol).unwrap();
        let js = logp::jensen_shannon_divergence(a, b, tol).unwrap();
        let tv = logp::total_variation(a, b, tol).unwrap();
        let bhatt_d = logp::bhattacharyya_distance(a, b, tol).unwrap();

        println!(
            "{:<14} {:>10.6} {:>10.6} {:>10.6} {:>10.6} {:>10.6} {:>10.6}",
            label, rao, hel, kl, js, tv, bhatt_d,
        );
    }

    println!();

    // --- Structural relationships ---
    // Rao distance is 2 * arccos(BC), where BC is the Bhattacharyya coefficient.
    // Hellinger distance is sqrt(1 - BC).
    // So both are monotone transformations of BC; verify consistency.
    println!("Consistency checks (Rao and Hellinger vs Bhattacharyya coefficient):");
    for &(label, a, b) in pairs {
        let bc = logp::bhattacharyya_coeff(a, b, tol).unwrap();
        let rao = rao_distance_categorical(a, b, tol).unwrap();
        let hel = hellinger(a, b, tol).unwrap();

        // Rao = 2 * arccos(BC)
        let rao_from_bc = 2.0 * bc.clamp(0.0, 1.0).acos();
        // Hellinger = sqrt(1 - BC)
        let hel_from_bc = (1.0 - bc).max(0.0).sqrt();

        let rao_ok = (rao - rao_from_bc).abs() < 1e-10;
        let hel_ok = (hel - hel_from_bc).abs() < 1e-10;

        println!("  {label:<14} BC={bc:.6}  Rao match={rao_ok}  Hellinger match={hel_ok}",);
        assert!(rao_ok, "Rao / BC mismatch for {label}");
        assert!(hel_ok, "Hellinger / BC mismatch for {label}");
    }

    println!();

    // --- Pinsker-type bound ---
    // TV(p,q) <= sqrt(KL(p||q) / 2)  (Pinsker's inequality).
    println!("Pinsker's inequality: TV <= sqrt(KL / 2)");
    for &(label, a, b) in pairs {
        let kl = logp::kl_divergence(a, b, tol).unwrap();
        let tv = logp::total_variation(a, b, tol).unwrap();
        let bound = (kl / 2.0).sqrt();
        println!(
            "  {label:<14} TV={tv:.6}  bound={bound:.6}  satisfied={}",
            tv <= bound + 1e-12
        );
        assert!(tv <= bound + 1e-12, "Pinsker violated for {label}");
    }

    println!();

    // --- Ordering ---
    // All divergences should agree: q_near is closer to p than q_far.
    let rao_near = rao_distance_categorical(&p, &q_near, tol).unwrap();
    let rao_far = rao_distance_categorical(&p, &q_far, tol).unwrap();
    let kl_near = logp::kl_divergence(&p, &q_near, tol).unwrap();
    let kl_far = logp::kl_divergence(&p, &q_far, tol).unwrap();
    let js_near = logp::jensen_shannon_divergence(&p, &q_near, tol).unwrap();
    let js_far = logp::jensen_shannon_divergence(&p, &q_far, tol).unwrap();

    println!("Ordering consistency (near < far):");
    println!(
        "  Rao:  {rao_near:.6} < {rao_far:.6} = {}",
        rao_near < rao_far
    );
    println!("  KL:   {kl_near:.6} < {kl_far:.6} = {}", kl_near < kl_far);
    println!("  JS:   {js_near:.6} < {js_far:.6} = {}", js_near < js_far);
    assert!(rao_near < rao_far);
    assert!(kl_near < kl_far);
    assert!(js_near < js_far);

    println!("\nAll checks passed.");
}
