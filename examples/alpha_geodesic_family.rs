//! The alpha-family of geodesics on the simplex, rendered in the terminal.
//!
//! Amari's alpha-geometry gives a one-parameter family of "straight lines"
//! between two distributions. The same pair of endpoints has a different
//! halfway point depending on which geometry you walk:
//!
//!   alpha = -1  mixture (m-) geodesic:    linear in probability    -> keeps both peaks
//!   alpha =  0  Fisher-Rao geodesic:      the Riemannian metric    -> in between
//!   alpha = +1  exponential (e-) geodesic: linear in log-prob       -> normalized geometric mean
//!
//! Walking two opposing humps to their t=0.5 midpoint makes the difference
//! visible: the mixture midpoint is the arithmetic average (bimodal), while the
//! exponential midpoint is the normalized geometric mean (unimodal, pulled to
//! the overlap). Fisher-Rao sits between them.
//!
//! Run: `cargo run --example alpha_geodesic_family`

use infogeom::alpha_geodesic;

fn main() {
    // Two opposing distributions on a 9-bin simplex: p peaks left, q peaks right.
    let p = bump(9, 2.0, 1.2);
    let q = bump(9, 6.0, 1.2);

    println!("=== Endpoints ===\n");
    println!("  p (peaks left):");
    bars(&p);
    println!("\n  q (peaks right):");
    bars(&q);

    println!("\n=== Midpoints (t = 0.5) under different geometries ===\n");
    for (alpha, name, note) in [
        (
            -1.0,
            "alpha = -1  (mixture)",
            "arithmetic average -- both peaks survive",
        ),
        (
            0.0,
            "alpha =  0  (Fisher-Rao)",
            "the Riemannian midpoint, between the two",
        ),
        (
            1.0,
            "alpha = +1  (exponential)",
            "normalized geometric mean -- pulled to the overlap",
        ),
    ] {
        let mid = alpha_geodesic(&p, &q, 0.5, alpha, 1e-12).expect("geodesic");
        println!("  {name}: {note}");
        bars(&mid);
        println!();
    }

    // Sweep alpha and report the heaviest bin of the midpoint. This is NOT
    // monotone in alpha: for alpha <= 0 the midpoint keeps the two ORIGINAL
    // modes (max bin sits at a side peak and the distribution is most spread
    // near alpha = 0), while for alpha > 0 mass concentrates at the central
    // overlap and the max bin climbs steeply. The two regimes peak in different
    // places, so the single max-bin number dips then rises.
    println!("=== Midpoint peakedness vs alpha (max bin mass) ===\n");
    for &alpha in &[-2.0, -1.0, 0.0, 1.0, 2.0] {
        let mid = alpha_geodesic(&p, &q, 0.5, alpha, 1e-12).expect("geodesic");
        let peak = mid.iter().copied().fold(0.0f64, f64::max);
        let n = (peak * 60.0).round() as usize;
        println!("  alpha={alpha:>4.1}  max={peak:.3}  {}", "#".repeat(n));
    }
}

// A discrete bump: bin i gets exp(-((i-mean)/std)^2 / 2), normalized.
fn bump(n: usize, mean: f64, std: f64) -> Vec<f64> {
    let mut v: Vec<f64> = (0..n)
        .map(|i| {
            let z = (i as f64 - mean) / std;
            (-0.5 * z * z).exp()
        })
        .collect();
    let s: f64 = v.iter().sum();
    for x in &mut v {
        *x /= s;
    }
    v
}

// Vertical bar chart of a probability vector (rows = height, columns = bins).
fn bars(p: &[f64]) {
    const ROWS: usize = 8;
    let max = p.iter().copied().fold(0.0f64, f64::max);
    let heights: Vec<usize> = p
        .iter()
        .map(|&x| {
            if max > 0.0 {
                (x / max * ROWS as f64).round() as usize
            } else {
                0
            }
        })
        .collect();
    for level in (1..=ROWS).rev() {
        let line: String = heights
            .iter()
            .map(|&h| if h >= level { " # " } else { "   " })
            .collect();
        println!("    {}", line.trim_end());
    }
    // bin index ruler
    let ruler: String = (0..p.len()).map(|i| format!("{i:^3}")).collect();
    println!("    {ruler}");
}
