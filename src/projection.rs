//! Fisher-orthogonal gradient projection for continual learning.
//!
//! Projects gradients onto the Fisher-orthogonal complement of stored past-task
//! gradients, preventing catastrophic forgetting while preserving current-task
//! performance.
//!
//! Given a current gradient `g`, past task gradients `G = [g_1, ..., g_k]`, and
//! a Fisher information matrix `F`, the projected gradient is:
//!
//! $$g_\text{proj} = g - G (G^\top F G)^{-1} G^\top F g$$
//!
//! This removes the component of `g` that lies in the Fisher-information-weighted
//! span of past gradients.
//!
//! Reference: arXiv 2601.12816, "Fisher-Orthogonal Projected Natural Gradient".

use crate::{Error, Result};

/// Project a gradient onto the Fisher-orthogonal complement of past gradients,
/// using a diagonal Fisher approximation.
///
/// This is the common case for large parameter spaces where storing the full
/// Fisher matrix is infeasible.
///
/// # Arguments
///
/// - `gradient`: current task gradient, length `d`
/// - `past_gradients`: `k` past task gradient vectors, each length `d`
/// - `fisher_diag`: diagonal of the Fisher information matrix, length `d`
///
/// # Returns
///
/// The projected gradient of length `d`, Fisher-orthogonal to all past gradients.
/// Returns the original gradient unchanged if `past_gradients` is empty.
pub fn fisher_orthogonal_project(
    gradient: &[f64],
    past_gradients: &[Vec<f64>],
    fisher_diag: &[f64],
) -> Result<Vec<f64>> {
    let d = gradient.len();
    validate_fisher_diag(fisher_diag, d)?;
    validate_past_gradients(past_gradients, d)?;

    if past_gradients.is_empty() {
        return Ok(gradient.to_vec());
    }

    let k = past_gradients.len();

    // Compute G^T F g  (k-vector)
    let gt_f_g: Vec<f64> = past_gradients
        .iter()
        .map(|gi| dot_diag_weighted(gi, fisher_diag, gradient))
        .collect();

    // Compute G^T F G  (k x k matrix, row-major)
    let mut gt_f_g_mat = vec![0.0; k * k];
    for i in 0..k {
        for j in i..k {
            let val = dot_diag_weighted(&past_gradients[i], fisher_diag, &past_gradients[j]);
            gt_f_g_mat[i * k + j] = val;
            gt_f_g_mat[j * k + i] = val;
        }
    }

    // Solve (G^T F G) x = G^T F g  via Cholesky
    let x = solve_symmetric_positive(&gt_f_g_mat, &gt_f_g, k)?;

    // g_proj = g - G x
    let mut result = gradient.to_vec();
    for (j, xj) in x.iter().enumerate() {
        for i in 0..d {
            result[i] -= past_gradients[j][i] * xj;
        }
    }

    Ok(result)
}

/// Project a gradient onto the Fisher-orthogonal complement of past gradients,
/// using a full Fisher information matrix.
///
/// For small parameter spaces where the full `d x d` Fisher matrix is available.
///
/// # Arguments
///
/// - `gradient`: current task gradient, length `d`
/// - `past_gradients`: `k` past task gradient vectors, each length `d`
/// - `fisher`: full Fisher information matrix, `d x d`, row-major
pub fn fisher_orthogonal_project_full(
    gradient: &[f64],
    past_gradients: &[Vec<f64>],
    fisher: &[f64],
) -> Result<Vec<f64>> {
    let d = gradient.len();

    if fisher.len() != d * d {
        return Err(Error::LengthMismatch {
            a_name: "fisher",
            a_len: fisher.len(),
            b_name: "gradient (d*d)",
            b_len: d * d,
        });
    }
    validate_past_gradients(past_gradients, d)?;

    if past_gradients.is_empty() {
        return Ok(gradient.to_vec());
    }

    let k = past_gradients.len();

    // F g  (d-vector)
    let f_g = matvec(fisher, gradient, d);

    // G^T F g  (k-vector)
    let gt_f_g: Vec<f64> = past_gradients.iter().map(|gi| dot(gi, &f_g)).collect();

    // F G_j for each past gradient (cache for reuse)
    let f_g_cols: Vec<Vec<f64>> = past_gradients
        .iter()
        .map(|gj| matvec(fisher, gj, d))
        .collect();

    // G^T F G  (k x k)
    let mut gt_f_g_mat = vec![0.0; k * k];
    for i in 0..k {
        for j in i..k {
            let val = dot(&past_gradients[i], &f_g_cols[j]);
            gt_f_g_mat[i * k + j] = val;
            gt_f_g_mat[j * k + i] = val;
        }
    }

    // Solve (G^T F G) x = G^T F g
    let x = solve_symmetric_positive(&gt_f_g_mat, &gt_f_g, k)?;

    // g_proj = g - G x
    let mut result = gradient.to_vec();
    for (j, xj) in x.iter().enumerate() {
        for i in 0..d {
            result[i] -= past_gradients[j][i] * xj;
        }
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Linear algebra helpers (small k, no external dep needed)
// ---------------------------------------------------------------------------

/// Dot product weighted by a diagonal matrix: sum_i a[i] * diag[i] * b[i].
fn dot_diag_weighted(a: &[f64], diag: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(diag.iter())
        .zip(b.iter())
        .map(|((&ai, &di), &bi)| ai * di * bi)
        .sum()
}

/// Standard dot product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&ai, &bi)| ai * bi).sum()
}

/// Matrix-vector product for a d x d row-major matrix.
fn matvec(mat: &[f64], v: &[f64], d: usize) -> Vec<f64> {
    (0..d)
        .map(|i| {
            let row = &mat[i * d..(i + 1) * d];
            dot(row, v)
        })
        .collect()
}

/// Solve A x = b for symmetric positive-definite A via Cholesky decomposition.
///
/// Falls back to regularization (A + eps * I) if the matrix is near-singular,
/// which happens when past gradients are nearly collinear in the Fisher metric.
fn solve_symmetric_positive(a: &[f64], b: &[f64], n: usize) -> Result<Vec<f64>> {
    // Cholesky: A = L L^T
    let mut l = vec![0.0; n * n];
    let eps = 1e-12;

    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i * n + k] * l[j * n + k];
            }
            if i == j {
                let diag = a[i * n + i] - sum + eps; // regularize
                if diag <= 0.0 {
                    // Matrix is not positive definite even with regularization.
                    // Return zero correction (no projection).
                    return Ok(vec![0.0; n]);
                }
                l[i * n + j] = diag.sqrt();
            } else {
                l[i * n + j] = (a[i * n + j] - sum) / l[j * n + j];
            }
        }
    }

    // Forward substitution: L y = b
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..i {
            sum += l[i * n + j] * y[j];
        }
        y[i] = (b[i] - sum) / l[i * n + i];
    }

    // Back substitution: L^T x = y
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..n {
            sum += l[j * n + i] * x[j];
        }
        x[i] = (y[i] - sum) / l[i * n + i];
    }

    Ok(x)
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

fn validate_fisher_diag(fisher_diag: &[f64], d: usize) -> Result<()> {
    if fisher_diag.len() != d {
        return Err(Error::LengthMismatch {
            a_name: "fisher_diag",
            a_len: fisher_diag.len(),
            b_name: "gradient",
            b_len: d,
        });
    }
    Ok(())
}

fn validate_past_gradients(past_gradients: &[Vec<f64>], d: usize) -> Result<()> {
    for (idx, g) in past_gradients.iter().enumerate() {
        if g.len() != d {
            return Err(Error::LengthMismatch {
                a_name: "past_gradient",
                a_len: g.len(),
                b_name: "gradient",
                b_len: d,
            });
        }
        let _ = idx; // suppress unused warning in non-debug
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    #[test]
    fn empty_past_returns_original() {
        let g = vec![1.0, 2.0, 3.0];
        let f = vec![1.0, 1.0, 1.0];
        let result = fisher_orthogonal_project(&g, &[], &f).unwrap();
        assert_eq!(result, g);
    }

    #[test]
    fn empty_past_returns_original_full() {
        let g = vec![1.0, 2.0];
        // 2x2 identity
        let f = vec![1.0, 0.0, 0.0, 1.0];
        let result = fisher_orthogonal_project_full(&g, &[], &f).unwrap();
        assert_eq!(result, g);
    }

    #[test]
    fn identity_fisher_reduces_to_standard_projection() {
        // With F = I, Fisher-orthogonal projection = standard orthogonal projection.
        let g = vec![3.0, 4.0];
        let past = vec![vec![1.0, 0.0]]; // past gradient along x-axis
        let f_diag = vec![1.0, 1.0];

        let proj = fisher_orthogonal_project(&g, &past, &f_diag).unwrap();

        // Standard projection onto orthogonal complement of [1,0]:
        // removes x-component, keeps y-component.
        assert!((proj[0]).abs() < EPS, "x should be 0, got {}", proj[0]);
        assert!(
            (proj[1] - 4.0).abs() < EPS,
            "y should be 4, got {}",
            proj[1]
        );
    }

    #[test]
    fn parallel_gradient_projects_to_zero() {
        let past = vec![vec![1.0, 2.0, 3.0]];
        let g = vec![2.0, 4.0, 6.0]; // parallel to past[0]
        let f_diag = vec![1.0, 1.0, 1.0];

        let proj = fisher_orthogonal_project(&g, &past, &f_diag).unwrap();

        for (i, &v) in proj.iter().enumerate() {
            assert!(v.abs() < EPS, "proj[{i}] = {v}, expected ~0");
        }
    }

    #[test]
    fn already_orthogonal_returns_unchanged() {
        // g = [0, 1], past = [[1, 0]], F = I
        // g is already orthogonal to past, so projection should be g itself.
        let g = vec![0.0, 1.0];
        let past = vec![vec![1.0, 0.0]];
        let f_diag = vec![1.0, 1.0];

        let proj = fisher_orthogonal_project(&g, &past, &f_diag).unwrap();

        assert!(
            (proj[0] - 0.0).abs() < EPS,
            "proj[0] = {}, expected 0",
            proj[0]
        );
        assert!(
            (proj[1] - 1.0).abs() < EPS,
            "proj[1] = {}, expected 1",
            proj[1]
        );
    }

    #[test]
    fn projected_is_fisher_orthogonal_to_past() {
        // Verify g_proj^T F g_past ~ 0 for each past gradient.
        let g = vec![1.0, 2.0, 3.0, 4.0];
        let past = vec![vec![1.0, 0.0, 1.0, 0.0], vec![0.0, 1.0, 0.0, 1.0]];
        let f_diag = vec![2.0, 3.0, 1.0, 4.0];

        let proj = fisher_orthogonal_project(&g, &past, &f_diag).unwrap();

        for (j, gj) in past.iter().enumerate() {
            let inner: f64 = proj
                .iter()
                .zip(f_diag.iter())
                .zip(gj.iter())
                .map(|((&pi, &fi), &gji)| pi * fi * gji)
                .sum();
            assert!(
                inner.abs() < EPS,
                "proj^T F g_past[{j}] = {inner}, expected ~0"
            );
        }
    }

    #[test]
    fn full_matrix_matches_diagonal_for_diagonal_fisher() {
        let g = vec![1.0, 2.0, 3.0];
        let past = vec![vec![1.0, 0.5, 0.0], vec![0.0, 0.5, 1.0]];
        let f_diag = vec![2.0, 3.0, 1.0];

        // Build full matrix from diagonal
        let mut f_full = vec![0.0; 9];
        for i in 0..3 {
            f_full[i * 3 + i] = f_diag[i];
        }

        let proj_diag = fisher_orthogonal_project(&g, &past, &f_diag).unwrap();
        let proj_full = fisher_orthogonal_project_full(&g, &past, &f_full).unwrap();

        for i in 0..3 {
            assert!(
                (proj_diag[i] - proj_full[i]).abs() < EPS,
                "diag[{i}]={} != full[{i}]={}",
                proj_diag[i],
                proj_full[i]
            );
        }
    }

    #[test]
    fn full_fisher_nondiagonal() {
        // Non-diagonal Fisher: F = [[2, 1], [1, 2]]
        let f = vec![2.0, 1.0, 1.0, 2.0];
        let g = vec![1.0, 1.0];
        let past = vec![vec![1.0, -1.0]];

        let proj = fisher_orthogonal_project_full(&g, &past, &f).unwrap();

        // Verify Fisher-orthogonality: proj^T F past[0] ~ 0
        let f_past: Vec<f64> = vec![
            f[0] * past[0][0] + f[1] * past[0][1],
            f[2] * past[0][0] + f[3] * past[0][1],
        ];
        let inner: f64 = proj[0] * f_past[0] + proj[1] * f_past[1];
        assert!(inner.abs() < EPS, "proj^T F past = {inner}, expected ~0");
    }

    #[test]
    fn length_mismatch_errors() {
        let g = vec![1.0, 2.0];
        let past = vec![vec![1.0, 2.0, 3.0]]; // wrong length
        let f_diag = vec![1.0, 1.0];

        assert!(fisher_orthogonal_project(&g, &past, &f_diag).is_err());

        let f_diag_wrong = vec![1.0, 1.0, 1.0]; // wrong length
        assert!(fisher_orthogonal_project(&g, &[], &f_diag_wrong).is_err());
    }

    #[test]
    fn multiple_past_gradients() {
        // 3D space, 2 past gradients spanning xy-plane, gradient along z.
        // With identity Fisher, projection should preserve the z-component.
        let g = vec![1.0, 1.0, 5.0];
        let past = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
        let f_diag = vec![1.0, 1.0, 1.0];

        let proj = fisher_orthogonal_project(&g, &past, &f_diag).unwrap();

        assert!(proj[0].abs() < EPS, "x should be 0, got {}", proj[0]);
        assert!(proj[1].abs() < EPS, "y should be 0, got {}", proj[1]);
        assert!(
            (proj[2] - 5.0).abs() < EPS,
            "z should be 5, got {}",
            proj[2]
        );
    }

    #[test]
    fn weighted_fisher_changes_projection() {
        // g = [1, 1], past = [[1, 1]].
        // With F = I: proj = 0 (g is parallel to past).
        // With F = [2, 1]: the "orthogonal complement" changes.
        let g = vec![1.0, 1.0];
        let past = vec![vec![1.0, 1.0]];

        let f_uniform = vec![1.0, 1.0];
        let proj_uniform = fisher_orthogonal_project(&g, &past, &f_uniform).unwrap();
        // Parallel: should be ~0
        for &v in &proj_uniform {
            assert!(v.abs() < EPS);
        }

        let f_weighted = vec![2.0, 1.0];
        let proj_weighted = fisher_orthogonal_project(&g, &past, &f_weighted).unwrap();
        // Still parallel in Fisher metric (g = past), so still ~0
        for &v in &proj_weighted {
            assert!(v.abs() < EPS);
        }

        // Now g not parallel to past in Fisher metric
        let g2 = vec![1.0, 0.0];
        let proj2 = fisher_orthogonal_project(&g2, &past, &f_weighted).unwrap();
        // Verify Fisher-orthogonality
        let inner: f64 = proj2
            .iter()
            .zip(f_weighted.iter())
            .zip(past[0].iter())
            .map(|((&pi, &fi), &gi)| pi * fi * gi)
            .sum();
        assert!(inner.abs() < EPS, "not Fisher-orthogonal: {inner}");
    }
}
