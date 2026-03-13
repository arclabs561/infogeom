//! Cross-crate example: optimal transport interpolation meets information geometry.
//!
//! Constructs a family of distributions via three interpolation strategies:
//! - Fisher-Rao geodesic (Riemannian geodesic on the simplex)
//! - Mixture (m-geodesic, alpha=-1)
//! - Displacement interpolation via Sinkhorn transport plan
//!
//! Then measures Fisher-Rao distance and Hellinger distance along each path from p0.
//!
//! Key insight: Fisher-Rao measures geodesic distance on the probability simplex
//! (intrinsic curvature), while Wasserstein/Sinkhorn measures the cost of mass
//! transport (ground-metric dependent). They can disagree on which paths are "short."

use infogeom::{fisher_rao_geodesic, hellinger, m_geodesic, rao_distance_categorical};
use ndarray::{Array1, Array2};

// ---------------------------------------------------------------------------
// Helpers: f32/f64 conversion between wass (f32) and infogeom (f64)
// ---------------------------------------------------------------------------

fn to_f32(v: &[f64]) -> Array1<f32> {
    Array1::from_iter(v.iter().map(|&x| x as f32))
}

fn to_f64_vec(v: &Array1<f32>) -> Vec<f64> {
    v.iter().map(|&x| f64::from(x)).collect()
}

/// Build a 1D Euclidean cost matrix for support points {0, 1, ..., n-1}.
fn cost_1d(n: usize) -> Array2<f32> {
    let mut c = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            c[[i, j]] = (i as f32 - j as f32).abs();
        }
    }
    c
}

/// Displacement interpolation via a transport plan T.
///
/// For each target bin j, accumulate the mass that T routes there, weighted by (1-t)
/// staying at origin i and t arriving at j. Concretely:
///
///   p_t[k] = sum_{i,j} T[i,j] * ((1-t)*delta(k,i) + t*delta(k,j))
///
/// This is the discrete McCann interpolation on a finite support.
fn displacement(p0: &[f64], plan: &Array2<f32>, t: f64, n: usize) -> Vec<f64> {
    let _ = p0; // p0 is encoded in the plan's row marginals
    let mut pt = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            let mass = f64::from(plan[[i, j]]);
            pt[i] += (1.0 - t) * mass;
            pt[j] += t * mass;
        }
    }
    // Renormalize to correct for small Sinkhorn residuals.
    let s: f64 = pt.iter().sum();
    if s > 0.0 {
        for x in &mut pt {
            *x /= s;
        }
    }
    pt
}

fn main() {
    let tol = 1e-12;
    let n = 5;

    // Two categorical distributions on {0,1,2,3,4}: p0 peaked left, p1 peaked right.
    let p0: Vec<f64> = vec![0.50, 0.30, 0.10, 0.05, 0.05];
    let p1: Vec<f64> = vec![0.05, 0.05, 0.10, 0.30, 0.50];

    // Sinkhorn transport plan between p0 and p1.
    let a = to_f32(&p0);
    let b = to_f32(&p1);
    let cost = cost_1d(n);
    let reg = 0.05_f32;
    let max_iter = 500;
    let (plan, ot_cost) = wass::sinkhorn_log(&a, &b, &cost, reg, max_iter);

    println!("Endpoint distributions:");
    println!("  p0 = {:?}", p0);
    println!("  p1 = {:?}", p1);
    println!("  Sinkhorn OT cost (reg={reg}): {ot_cost:.6}",);
    println!();

    // --- Table header ---
    println!(
        "{:>4}  {:>12} {:>12}  {:>12} {:>12}  {:>12} {:>12}",
        "t", "FR(fr-geo)", "FR(mix)", "FR(disp)", "H(fr-geo)", "H(mix)", "H(disp)"
    );
    println!("{}", "-".repeat(92));

    let steps: Vec<f64> = (0..=10).map(|i| i as f64 / 10.0).collect();

    for &t in &steps {
        let pfr = fisher_rao_geodesic(&p0, &p1, t, tol).expect("fisher-rao geodesic");
        let pm = m_geodesic(&p0, &p1, t, tol).expect("m-geodesic");
        let pd = displacement(&p0, &plan, t, n);

        // Fisher-Rao distance from p0 to interpolated distribution.
        let fr_geo = rao_distance_categorical(&p0, &pfr, tol).expect("rao (fr-geo)");
        let fr_mix = rao_distance_categorical(&p0, &pm, tol).expect("rao (mix)");
        let fr_disp = rao_distance_categorical(&p0, &pd, tol).expect("rao (disp)");

        // Hellinger distance from p0.
        let h_geo = hellinger(&p0, &pfr, tol).expect("hellinger (fr-geo)");
        let h_mix = hellinger(&p0, &pm, tol).expect("hellinger (mix)");
        let h_disp = hellinger(&p0, &pd, tol).expect("hellinger (disp)");

        println!(
            "{t:>4.1}  {fr_geo:>12.6} {fr_mix:>12.6}  {fr_disp:>12.6} {h_geo:>12.6}  {h_mix:>12.6} {h_disp:>12.6}",
        );
    }

    println!();

    // --- Verification: endpoint sanity ---
    let pm0 = m_geodesic(&p0, &p1, 0.0, tol).expect("m-geo t=0");
    let pd0 = displacement(&p0, &plan, 0.0, n);
    assert!(
        rao_distance_categorical(&p0, &pm0, tol).expect("rao t=0 mix") < 1e-10,
        "mixture at t=0 should equal p0"
    );
    assert!(
        rao_distance_categorical(&p0, &pd0, tol).expect("rao t=0 disp") < 1e-6,
        "displacement at t=0 should be near p0"
    );

    let fr_full = rao_distance_categorical(&p0, &p1, tol).expect("rao full");
    let pm1 = m_geodesic(&p0, &p1, 1.0, tol).expect("m-geo t=1");
    let fr_mix_1 = rao_distance_categorical(&p0, &pm1, tol).expect("rao t=1");
    assert!(
        (fr_full - fr_mix_1).abs() < 1e-10,
        "mixture at t=1 should equal p1"
    );

    // --- Metric comparison ---
    // Compute cumulative arc lengths along each path for Fisher-Rao.
    let mut arc_geo = 0.0_f64;
    let mut arc_mix = 0.0_f64;
    let mut arc_disp = 0.0_f64;
    for i in 1..steps.len() {
        let prev_fr = fisher_rao_geodesic(&p0, &p1, steps[i - 1], tol).expect("arc fr prev");
        let curr_fr = fisher_rao_geodesic(&p0, &p1, steps[i], tol).expect("arc fr curr");
        arc_geo += rao_distance_categorical(&prev_fr, &curr_fr, tol).expect("arc fr");

        let prev_m = m_geodesic(&p0, &p1, steps[i - 1], tol).expect("arc mix prev");
        let curr_m = m_geodesic(&p0, &p1, steps[i], tol).expect("arc mix curr");
        arc_mix += rao_distance_categorical(&prev_m, &curr_m, tol).expect("arc mix");

        let prev_d = displacement(&p0, &plan, steps[i - 1], n);
        let curr_d = displacement(&p0, &plan, steps[i], n);
        arc_disp += rao_distance_categorical(&prev_d, &curr_d, tol).expect("arc disp");
    }

    println!("Cumulative Fisher-Rao arc length (sum of consecutive distances):");
    println!("  FR geodesic path:   {arc_geo:.6}");
    println!("  Mixture path:       {arc_mix:.6}");
    println!("  Displacement path:  {arc_disp:.6}");
    println!("  Direct p0->p1:      {fr_full:.6}");
    println!();
    println!("The direct distance is a lower bound (triangle inequality).");
    println!("The FR geodesic path should have arc length closest to the direct distance.");
    println!("Mixture interpolation is NOT the Fisher-Rao geodesic,");
    println!("so its arc length exceeds the direct distance.");

    // Sinkhorn divergence at endpoints for comparison.
    let s_full = wass::sinkhorn_divergence_same_support(&a, &b, &cost, reg, max_iter, 1e-6)
        .expect("sinkhorn div full");
    println!();
    println!("Sinkhorn divergence p0 vs p1: {s_full:.6}");

    // Verify triangle inequality for Fisher-Rao arc lengths.
    assert!(
        arc_geo + 1e-8 >= fr_full,
        "triangle inequality violated for FR geodesic path"
    );
    assert!(
        arc_mix + 1e-8 >= fr_full,
        "triangle inequality violated for mixture path"
    );
    assert!(
        arc_disp + 1e-8 >= fr_full,
        "triangle inequality violated for displacement path"
    );

    // --- f64/f32 round-trip sanity ---
    let roundtrip: Vec<f64> = to_f64_vec(&to_f32(&p0));
    let max_err: f64 = p0
        .iter()
        .zip(roundtrip.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        max_err < 1e-6,
        "f64->f32->f64 round-trip error too large: {max_err}"
    );

    println!();
    println!("All checks passed.");
}
