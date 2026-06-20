//! Detecting retrieval distribution shift on the probability simplex.
//!
//! A ranker often produces a distribution over result families: API docs, examples,
//! design notes, issue threads, and unrelated hits. When the distribution moves,
//! Euclidean distance treats all coordinates as flat dimensions. Fisher-Rao distance
//! measures the move on the simplex itself.
//!
//! Run:
//!
//! ```sh
//! cargo run --example retrieval_distribution_shift
//! ```

use infogeom::{
    fisher_rao_geodesic, hellinger, m_geodesic, natural_gradient, rao_distance_categorical,
};

const LABELS: [&str; 5] = ["api", "examples", "design", "issues", "noise"];

fn print_distribution(label: &str, p: &[f64]) {
    println!("{label}");
    for (name, value) in LABELS.iter().zip(p.iter()) {
        let bars = (value * 40.0).round() as usize;
        println!("  {name:>8}: {value:>5.3}  {}", "#".repeat(bars));
    }
}

fn main() {
    let tol = 1e-12;

    // A query before and after adding a new document family to the index.
    // The post-update distribution has moved from API/examples toward design notes
    // and issue threads, with a little more noise.
    let before = [0.44, 0.30, 0.16, 0.07, 0.03];
    let after = [0.16, 0.18, 0.38, 0.20, 0.08];

    print_distribution("before", &before);
    println!();
    print_distribution("after", &after);
    println!();

    let rao = rao_distance_categorical(&before, &after, tol).expect("rao distance");
    let h = hellinger(&before, &after, tol).expect("hellinger");
    println!("Shift metrics");
    println!("  Fisher-Rao distance: {rao:.6} radians");
    println!("  Hellinger distance:  {h:.6}");
    println!();

    let fr_mid = fisher_rao_geodesic(&before, &after, 0.5, tol).expect("fr midpoint");
    let mix_mid = m_geodesic(&before, &after, 0.5, tol).expect("mixture midpoint");

    println!("Midpoint comparison");
    println!("{:>8}  {:>10}  {:>10}", "family", "Fisher-Rao", "mixture");
    println!("{}", "-".repeat(34));
    for ((name, fr), mix) in LABELS.iter().zip(fr_mid.iter()).zip(mix_mid.iter()) {
        println!("{name:>8}  {fr:>10.4}  {mix:>10.4}");
    }
    println!();

    let d_left = rao_distance_categorical(&before, &fr_mid, tol).expect("left half");
    let d_right = rao_distance_categorical(&fr_mid, &after, tol).expect("right half");
    println!("Fisher-Rao midpoint check");
    println!("  d(before, midpoint): {d_left:.6}");
    println!("  d(midpoint, after):  {d_right:.6}");
    println!("  d(before, after):    {rao:.6}");
    assert!((d_left - d_right).abs() < 1e-10);
    assert!((d_left + d_right - rao).abs() < 1e-10);
    println!();

    // Suppose human feedback says examples and design notes are useful, noise and
    // issue threads are not. Natural gradient scales the Euclidean feedback by the
    // current distribution, so tiny coordinates do not dominate the update.
    let feedback_grad = [-0.15, 0.45, 0.35, -0.20, -0.45];
    let ng = natural_gradient(&after, &feedback_grad).expect("natural gradient");

    println!("Feedback step at the shifted distribution");
    println!("{:>8}  {:>10}  {:>10}", "family", "euclidean", "natural");
    println!("{}", "-".repeat(34));
    for ((name, g), nat) in LABELS.iter().zip(feedback_grad.iter()).zip(ng.iter()) {
        println!("{name:>8}  {g:>10.4}  {nat:>10.4}");
    }

    println!();
    println!("The natural step is small on the rare noise bucket even though the raw");
    println!("feedback there is large. That is the inverse Fisher metric at work.");

    assert!(rao > 0.0);
    assert!(h > 0.0);
    assert!(ng[1] > 0.0);
    assert!(ng[4] < 0.0);
}
