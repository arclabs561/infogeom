//! Fisher-orthogonal projection must hold at small gradient and Fisher scales.

use infogeom::*;

fn fdot(a: &[f64], f: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(f).zip(b).map(|((x, w), y)| x * w * y).sum()
}

/// Docs: "The projected gradient ... Fisher-orthogonal to all past
/// gradients." The Cholesky adds an absolute 1e-12 ridge, so when
/// G^T F G is ~1e-12 or smaller (small gradients or a small empirical Fisher
/// diagonal) most of the past-gradient component survives.
#[test]
fn projection_is_fisher_orthogonal_at_small_scale() {
    let g = [1e-6, 2e-6, -1e-6];
    let past = vec![vec![1e-6, 0.0, 0.0]];
    let f = [1.0, 1.0, 1.0];
    let out = fisher_orthogonal_project(&g, &past, &f).unwrap();
    let ip = fdot(&out, &f, &past[0]);
    let scale = fdot(&g, &f, &past[0]).abs();
    assert!(
        ip.abs() <= 1e-6 * scale,
        "<proj, past>_F = {ip:e} (unprojected {scale:e}); out={out:?}"
    );
}

#[test]
fn projection_is_fisher_orthogonal_with_small_fisher_diag() {
    // EWC-style diagonal Fisher = mean squared gradient ~ 1e-8.
    let g = [1e-3, 1e-3];
    let past = vec![vec![1e-3, 0.0]];
    let f = [1e-8, 1e-8];
    let out = fisher_orthogonal_project(&g, &past, &f).unwrap();
    assert!(
        out[0].abs() < 1e-9,
        "component along past survived: {out:?}"
    );
}

/// Same scale problem in the full-Fisher variant.
#[test]
fn projection_full_is_fisher_orthogonal_at_small_scale() {
    let g = [1e-6, 1e-6];
    let past = vec![vec![1e-6, 0.0]];
    let f = [1.0, 0.0, 0.0, 1.0];
    let out = fisher_orthogonal_project_full(&g, &past, &f).unwrap();
    assert!(out[0].abs() < 1e-12, "out={out:?}");
}
