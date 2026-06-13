//! `infogeom`: information geometry on the probability simplex.
//!
//! Provides geodesic interpolation, natural gradient, and Fisher information
//! for categorical distributions. Builds on [`logp`] for divergence primitives.
//!
//! ## Geodesics
//!
//! Three geodesic families on the simplex, corresponding to the alpha-connection
//! framework of Amari & Nagaoka (2000, Ch. 3):
//!
//! - [`fisher_rao_geodesic`]: the unique Riemannian (alpha=0) geodesic via the sphere embedding
//! - [`m_geodesic`]: the mixture (alpha=-1) geodesic -- linear interpolation in probability space
//! - [`e_geodesic`]: the exponential (alpha=+1) geodesic -- linear interpolation in log space
//!
//! ## Distances
//!
//! - [`rao_distance_categorical`]: Fisher-Rao distance (radians) via the Bhattacharyya coefficient
//! - Hellinger distance: re-exported from `logp` as [`hellinger`]
//!
//! ## Fisher information
//!
//! - [`fisher_information_diagonal`]: diagonal of the Fisher information matrix for categoricals
//! - [`natural_gradient`]: multiply Euclidean gradient by the inverse Fisher metric
//!
//! ## Fisher-orthogonal projection
//!
//! - [`fisher_orthogonal_project`]: project gradient onto Fisher-orthogonal complement
//!   of past task gradients (diagonal Fisher approximation)
//! - [`fisher_orthogonal_project_full`]: same, with full Fisher matrix
//!
//! ## Quick example
//!
//! ```rust
//! use infogeom::{fisher_rao_geodesic, rao_distance_categorical};
//!
//! let p = [0.7, 0.2, 0.1];
//! let q = [0.1, 0.2, 0.7];
//!
//! let mid = fisher_rao_geodesic(&p, &q, 0.5, 1e-12).unwrap();
//! let d_full = rao_distance_categorical(&p, &q, 1e-12).unwrap();
//! let d_half = rao_distance_categorical(&p, &mid, 1e-12).unwrap();
//!
//! // Midpoint is equidistant from both endpoints.
//! assert!((d_half - d_full / 2.0).abs() < 1e-10);
//! ```
//!
//! # Related crates
//! - [`qig`]: Quantum generalization of Fisher-Rao geometry to density matrices.
//! - [`logp`]: Provides the underlying divergences (KL, Hellinger) that infogeom's geodesics are built from.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

#[cfg(feature = "manifold")]
mod manifold;
#[cfg(feature = "manifold")]
pub use manifold::FisherRaoSimplex;

mod projection;
pub use projection::{fisher_orthogonal_project, fisher_orthogonal_project_full};

use thiserror::Error;

/// Re-export Hellinger distance from `logp`.
pub use logp::hellinger;

/// Errors produced by `infogeom` operations.
#[derive(Debug, Error)]
pub enum Error {
    /// Propagated from `logp` (simplex validation, etc.).
    #[error(transparent)]
    Logp(#[from] logp::Error),

    /// Interpolation parameter `t` is outside [0, 1].
    #[error("interpolation parameter t={t} outside [0, 1]")]
    InvalidT {
        /// The out-of-range value.
        t: f64,
    },

    /// Two inputs have incompatible lengths.
    #[error("length mismatch: {a_name} has {a_len} elements, {b_name} has {b_len}")]
    LengthMismatch {
        /// Name of the first argument.
        a_name: &'static str,
        /// Length of the first argument.
        a_len: usize,
        /// Name of the second argument.
        b_name: &'static str,
        /// Length of the second argument.
        b_len: usize,
    },

    /// A zero entry was found where strictly positive values are required.
    #[error(
        "zero entry at index {idx}: log-space operations require strictly positive distributions"
    )]
    ZeroEntry {
        /// Index of the zero (or negative) entry.
        idx: usize,
    },
}

/// Crate-level result type.
pub type Result<T> = core::result::Result<T, Error>;

// ---------------------------------------------------------------------------
// Validation helpers
// ---------------------------------------------------------------------------

/// Validate that `t` is in [0, 1].
fn validate_t(t: f64) -> Result<()> {
    if !(0.0..=1.0).contains(&t) {
        return Err(Error::InvalidT { t });
    }
    Ok(())
}

/// Validate that two slices have the same length, returning a descriptive error.
fn validate_same_len(
    a: &[f64],
    b: &[f64],
    a_name: &'static str,
    b_name: &'static str,
) -> Result<()> {
    if a.len() != b.len() {
        return Err(Error::LengthMismatch {
            a_name,
            a_len: a.len(),
            b_name,
            b_len: b.len(),
        });
    }
    Ok(())
}

/// Validate that all entries are strictly positive, returning `ZeroEntry` on the first violation.
fn validate_strictly_positive(p: &[f64]) -> Result<()> {
    for (i, &v) in p.iter().enumerate() {
        if v <= 0.0 {
            return Err(Error::ZeroEntry { idx: i });
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Distances
// ---------------------------------------------------------------------------

/// Fisher-Rao (Rao) distance between categorical distributions `p` and `q` (in radians).
///
/// For categorical distributions, the Fisher-Rao manifold embeds isometrically onto a
/// sphere via the map p -> 2*sqrt(p). Under this embedding the geodesic distance is:
///
///   d_FR(p, q) = 2 * arccos(BC(p, q))
///
/// where BC is the Bhattacharyya coefficient.
///
/// # Snap heuristic
///
/// When `BC(p,q)` is within `10 * tol` of 1.0, we snap it to exactly 1.0 so that
/// `d(p, p) == 0` despite floating-point error. The factor 10x accounts for accumulated
/// rounding in the BC summation: each of the n sqrt-multiply-accumulate steps can
/// introduce up to ~eps_machine error, so the total deviation from the true BC is
/// bounded roughly by n * eps. Using 10 * tol provides margin for dimensions up to
/// approximately 10 / eps_machine (about 10^17 for f64 with tol=1e-12).
///
/// Reference: Amari & Nagaoka (2000), Ch. 2, Thm. 2.2 (Fisher metric on statistical
/// models).
pub fn rao_distance_categorical(p: &[f64], q: &[f64], tol: f64) -> Result<f64> {
    let mut bc = logp::bhattacharyya_coeff(p, q, tol)?;

    bc = bc.clamp(0.0, 1.0);
    // Snap BC to 1.0 when within 10*tol of 1.0. See doc comment above for rationale.
    if (1.0 - bc).abs() <= 10.0 * tol {
        bc = 1.0;
    }
    Ok(2.0 * bc.acos())
}

// ---------------------------------------------------------------------------
// Geodesics
// ---------------------------------------------------------------------------

/// Fisher-Rao geodesic (slerp on the simplex via the sphere embedding).
///
/// The categorical simplex embeds isometrically into the positive orthant of the unit
/// sphere via theta_i = sqrt(p_i). The geodesic on the sphere is the great-circle
/// (slerp), pulled back to the simplex by squaring:
///
///   theta_i(t) = (sin((1-t)*omega) * sqrt(p_i) + sin(t*omega) * sqrt(q_i)) / sin(omega)
///   gamma_i(t) = theta_i(t)^2
///
/// where omega = arccos(BC(p, q)).
///
/// This is the alpha=0 geodesic in Amari's alpha-connection framework.
///
/// Reference: Amari & Nagaoka (2000), Ch. 2 (sphere embedding of the categorical
/// Fisher-Rao manifold).
pub fn fisher_rao_geodesic(p: &[f64], q: &[f64], t: f64, tol: f64) -> Result<Vec<f64>> {
    validate_t(t)?;
    logp::validate_simplex(p, tol)?;
    logp::validate_simplex(q, tol)?;
    validate_same_len(p, q, "p", "q")?;

    let bc = logp::bhattacharyya_coeff(p, q, tol)?.clamp(0.0, 1.0);
    let omega = bc.acos();

    // Degenerate case: p and q are (nearly) identical.
    if omega < tol {
        return Ok(p.to_vec());
    }

    let sin_omega = omega.sin();
    let s0 = ((1.0 - t) * omega).sin();
    let s1 = (t * omega).sin();

    let result: Vec<f64> = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| {
            let theta = (s0 * pi.sqrt() + s1 * qi.sqrt()) / sin_omega;
            theta * theta
        })
        .collect();

    Ok(result)
}

/// Mixture geodesic (alpha = -1): linear interpolation in probability space.
///
///   gamma_i(t) = (1 - t) * p_i + t * q_i
///
/// This always stays on the simplex (convex combination of simplices).
///
/// Reference: Amari & Nagaoka (2000), Ch. 3, Sec. 3.4 (m-connection).
pub fn m_geodesic(p: &[f64], q: &[f64], t: f64, tol: f64) -> Result<Vec<f64>> {
    validate_t(t)?;
    logp::validate_simplex(p, tol)?;
    logp::validate_simplex(q, tol)?;
    validate_same_len(p, q, "p", "q")?;

    let result: Vec<f64> = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (1.0 - t) * pi + t * qi)
        .collect();

    Ok(result)
}

/// Exponential geodesic (alpha = +1): linear interpolation in log space.
///
///   log_unnorm_i = (1 - t) * ln(p_i) + t * ln(q_i)
///   gamma_i(t) = exp(log_unnorm_i) / sum_j exp(log_unnorm_j)
///
/// Requires all entries of p and q to be strictly positive (returns [`Error::ZeroEntry`]
/// otherwise). Uses log-sum-exp for numerical stability.
///
/// Reference: Amari & Nagaoka (2000), Ch. 3, Sec. 3.4 (e-connection).
pub fn e_geodesic(p: &[f64], q: &[f64], t: f64, tol: f64) -> Result<Vec<f64>> {
    validate_t(t)?;
    logp::validate_simplex(p, tol)?;
    logp::validate_simplex(q, tol)?;
    validate_same_len(p, q, "p", "q")?;
    validate_strictly_positive(p)?;
    validate_strictly_positive(q)?;

    // Compute unnormalized log-probabilities.
    let log_unnorm: Vec<f64> = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (1.0 - t) * pi.ln() + t * qi.ln())
        .collect();

    // Log-sum-exp for numerical stability.
    let lse = logp::log_sum_exp(&log_unnorm);

    let result: Vec<f64> = log_unnorm.iter().map(|&v| (v - lse).exp()).collect();

    Ok(result)
}

/// Alpha-geodesic on the simplex (Amari alpha-connection), unifying the
/// mixture, exponential, and Fisher-Rao families in one parametric surface.
///
/// Affine parameterization in the alpha-embedding `f_alpha(p)` proportional to
/// `p^((1 - alpha) / 2)`: interpolate the embeddings linearly, then renormalize
/// onto the simplex.
///
/// - `alpha = -1.0` reduces exactly to [`m_geodesic`] (linear in `p`).
/// - `alpha = 1.0` reduces exactly to [`e_geodesic`] (linear in `log p`).
/// - `alpha = 0.0` traces the same curve as [`fisher_rao_geodesic`] (the
///   sqrt-embedding great circle); the curve is identical, but the parameter is
///   the affine-alpha one rather than the arc-length (slerp) one, so the spacing
///   of intermediate points differs.
///
/// `alpha > 1.0` uses a negative embedding exponent and so requires strictly
/// positive entries (same zero-entry error as the exponential geodesic);
/// `alpha <= 1.0` allows zeros.
///
/// Reference: Amari & Nagaoka (2000), Ch. 2-3 (alpha-connections).
pub fn alpha_geodesic(p: &[f64], q: &[f64], t: f64, alpha: f64, tol: f64) -> Result<Vec<f64>> {
    validate_t(t)?;
    logp::validate_simplex(p, tol)?;
    logp::validate_simplex(q, tol)?;
    validate_same_len(p, q, "p", "q")?;

    // alpha = +1 is the log-embedding limit; delegate to the stable e-geodesic.
    if (alpha - 1.0).abs() < tol {
        return e_geodesic(p, q, t, tol);
    }

    let m = (1.0 - alpha) / 2.0; // alpha-embedding exponent
    if m < 0.0 {
        validate_strictly_positive(p)?;
        validate_strictly_positive(q)?;
    }

    let unnorm: Vec<f64> = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| {
            let mixed = (1.0 - t) * pi.powf(m) + t * qi.powf(m);
            mixed.powf(1.0 / m)
        })
        .collect();

    let sum: f64 = unnorm.iter().sum();
    Ok(unnorm.iter().map(|&v| v / sum).collect())
}

// ---------------------------------------------------------------------------
// Fisher information and natural gradient
// ---------------------------------------------------------------------------

/// Fisher information matrix diagonal for a categorical distribution.
///
/// For the categorical family, the Fisher information matrix is diagonal:
/// F_ij = delta_ij / p_i. This function returns the diagonal [1/p_1, ..., 1/p_n].
///
/// Returns [`Error::ZeroEntry`] if any p_i <= 0 (Fisher information is undefined at the
/// simplex boundary).
///
/// Reference: Amari & Nagaoka (2000), Ch. 2, Example 2.1.
pub fn fisher_information_diagonal(p: &[f64], tol: f64) -> Result<Vec<f64>> {
    logp::validate_simplex(p, tol)?;
    validate_strictly_positive(p)?;
    Ok(p.iter().map(|&pi| 1.0 / pi).collect())
}

/// Natural gradient on the categorical manifold.
///
/// The Fisher information matrix for categorical distributions is diag(1/p_i).
/// Its inverse is diag(p_i), so the natural gradient is the elementwise product:
///
///   (nabla_nat f)_i = p_i * g_i
///
/// where g is the Euclidean gradient.
///
/// Reference: Amari (1998), "Natural Gradient Works Efficiently in Learning",
/// Neural Computation 10(2), pp. 251-276.
pub fn natural_gradient(p: &[f64], euclidean_grad: &[f64]) -> Result<Vec<f64>> {
    validate_same_len(p, euclidean_grad, "p", "euclidean_grad")?;
    Ok(p.iter()
        .zip(euclidean_grad.iter())
        .map(|(&pi, &gi)| pi * gi)
        .collect())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Generate a random probability vector on the simplex (may have zeros).
    fn simplex_vec(len: usize) -> impl Strategy<Value = Vec<f64>> {
        prop::collection::vec(0.0f64..10.0, len).prop_map(|mut v| {
            let s: f64 = v.iter().sum();
            if s == 0.0 {
                v[0] = 1.0;
                return v;
            }
            for x in v.iter_mut() {
                *x /= s;
            }
            v
        })
    }

    /// Generate a strictly positive probability vector (no zeros).
    fn positive_simplex_vec(len: usize) -> impl Strategy<Value = Vec<f64>> {
        prop::collection::vec(0.01f64..10.0, len).prop_map(|mut v| {
            let s: f64 = v.iter().sum();
            for x in v.iter_mut() {
                *x /= s;
            }
            v
        })
    }

    const TOL: f64 = 1e-9;

    // -- Existing distance tests --

    proptest! {
        #[test]
        fn rao_is_symmetric(p in simplex_vec(8), q in simplex_vec(8)) {
            let d1 = rao_distance_categorical(&p, &q, TOL).unwrap();
            let d2 = rao_distance_categorical(&q, &p, TOL).unwrap();
            prop_assert!((d1 - d2).abs() < 1e-12);
        }

        #[test]
        fn rao_self_is_zero(p in simplex_vec(12)) {
            let d = rao_distance_categorical(&p, &p, TOL).unwrap();
            prop_assert!(d.abs() < 1e-12);
        }

        #[test]
        fn rao_is_bounded(p in simplex_vec(10), q in simplex_vec(10)) {
            let d = rao_distance_categorical(&p, &q, TOL).unwrap();
            prop_assert!(d >= -1e-12);
            prop_assert!(d <= core::f64::consts::PI + 1e-12);
        }

        #[test]
        fn hellinger_is_symmetric(p in simplex_vec(8), q in simplex_vec(8)) {
            let h1 = hellinger(&p, &q, TOL).unwrap();
            let h2 = hellinger(&q, &p, TOL).unwrap();
            prop_assert!((h1 - h2).abs() < 1e-12);
        }

        #[test]
        fn hellinger_is_bounded(p in simplex_vec(8), q in simplex_vec(8)) {
            let h = hellinger(&p, &q, TOL).unwrap();
            prop_assert!(h >= -1e-12);
            prop_assert!(h <= 1.0 + 1e-12);
        }

        #[test]
        fn rao_matches_bhattacharyya_formula(p in simplex_vec(10), q in simplex_vec(10)) {
            let rao = rao_distance_categorical(&p, &q, TOL).unwrap();
            let mut bc = logp::bhattacharyya_coeff(&p, &q, TOL).unwrap();
            bc = bc.clamp(0.0, 1.0);
            if (1.0 - bc).abs() <= 10.0 * TOL {
                bc = 1.0;
            }
            let expected = 2.0 * bc.acos();
            prop_assert!((rao - expected).abs() < 1e-12, "rao={rao} expected={expected} bc={bc}");
        }
    }

    // -- New tests --

    proptest! {
        /// Triangle inequality: d(p,r) <= d(p,q) + d(q,r).
        #[test]
        fn rao_triangle_inequality(
            p in simplex_vec(5),
            q in simplex_vec(5),
            r in simplex_vec(5),
        ) {
            let dpq = rao_distance_categorical(&p, &q, TOL).unwrap();
            let dqr = rao_distance_categorical(&q, &r, TOL).unwrap();
            let dpr = rao_distance_categorical(&p, &r, TOL).unwrap();
            prop_assert!(dpr <= dpq + dqr + 1e-10,
                "triangle inequality violated: d(p,r)={dpr} > d(p,q)+d(q,r)={}", dpq + dqr);
        }
    }

    #[test]
    fn rao_disjoint_equals_pi() {
        // Disjoint distributions on the simplex: maximum distance = pi.
        let p = vec![1.0, 0.0, 0.0, 0.0];
        let q = vec![0.0, 1.0, 0.0, 0.0];
        let d = rao_distance_categorical(&p, &q, 1e-12).unwrap();
        assert!(
            (d - core::f64::consts::PI).abs() < 1e-10,
            "expected pi, got {d}"
        );
    }

    proptest! {
        /// hellinger(p, p) == 0 for any valid simplex vector.
        /// Tolerance is 1e-8 rather than 1e-12 because the BC summation
        /// accumulates rounding error proportional to dimension.
        #[test]
        fn hellinger_self_distance_zero(p in simplex_vec(8)) {
            let h = hellinger(&p, &p, TOL).unwrap();
            prop_assert!(h.abs() < 1e-7, "hellinger(p,p) = {h}, expected 0");
        }
    }

    // -- Geodesic endpoint tests --

    proptest! {
        /// Fisher-Rao geodesic at t=0 returns p, at t=1 returns q.
        #[test]
        fn fisher_rao_geodesic_endpoints(
            p in positive_simplex_vec(6),
            q in positive_simplex_vec(6),
        ) {
            let g0 = fisher_rao_geodesic(&p, &q, 0.0, TOL).unwrap();
            let g1 = fisher_rao_geodesic(&p, &q, 1.0, TOL).unwrap();

            for i in 0..p.len() {
                prop_assert!((g0[i] - p[i]).abs() < 1e-10,
                    "endpoint t=0: g0[{i}]={} != p[{i}]={}", g0[i], p[i]);
                prop_assert!((g1[i] - q[i]).abs() < 1e-10,
                    "endpoint t=1: g1[{i}]={} != q[{i}]={}", g1[i], q[i]);
            }
        }
    }

    proptest! {
        /// The midpoint of the Fisher-Rao geodesic should satisfy:
        /// d(p, mid) + d(mid, q) ~= d(p, q).
        #[test]
        fn fisher_rao_geodesic_midpoint_minimizes_sum_distances(
            p in positive_simplex_vec(5),
            q in positive_simplex_vec(5),
        ) {
            let mid = fisher_rao_geodesic(&p, &q, 0.5, TOL).unwrap();
            let d_full = rao_distance_categorical(&p, &q, TOL).unwrap();
            let d_pm = rao_distance_categorical(&p, &mid, TOL).unwrap();
            let d_mq = rao_distance_categorical(&mid, &q, TOL).unwrap();

            prop_assert!((d_pm + d_mq - d_full).abs() < 1e-8,
                "d(p,mid)+d(mid,q)={} != d(p,q)={}", d_pm + d_mq, d_full);
        }
    }

    proptest! {
        /// Mixture geodesic at t=0 returns p, at t=1 returns q.
        #[test]
        fn m_geodesic_endpoints(
            p in simplex_vec(6),
            q in simplex_vec(6),
        ) {
            let g0 = m_geodesic(&p, &q, 0.0, TOL).unwrap();
            let g1 = m_geodesic(&p, &q, 1.0, TOL).unwrap();

            for i in 0..p.len() {
                prop_assert!((g0[i] - p[i]).abs() < 1e-14,
                    "m-geodesic t=0: g0[{i}]={} != p[{i}]={}", g0[i], p[i]);
                prop_assert!((g1[i] - q[i]).abs() < 1e-14,
                    "m-geodesic t=1: g1[{i}]={} != q[{i}]={}", g1[i], q[i]);
            }
        }
    }

    proptest! {
        /// Exponential geodesic at t=0 returns p, at t=1 returns q.
        /// Requires strictly positive distributions.
        #[test]
        fn e_geodesic_endpoints(
            p in positive_simplex_vec(6),
            q in positive_simplex_vec(6),
        ) {
            let g0 = e_geodesic(&p, &q, 0.0, TOL).unwrap();
            let g1 = e_geodesic(&p, &q, 1.0, TOL).unwrap();

            for i in 0..p.len() {
                prop_assert!((g0[i] - p[i]).abs() < 1e-10,
                    "e-geodesic t=0: g0[{i}]={} != p[{i}]={}", g0[i], p[i]);
                prop_assert!((g1[i] - q[i]).abs() < 1e-10,
                    "e-geodesic t=1: g1[{i}]={} != q[{i}]={}", g1[i], q[i]);
            }
        }
    }

    proptest! {
        /// Along the e-geodesic from p to q, KL(p || gamma(t)) is monotonically
        /// non-decreasing in t. This follows from the convexity of KL in its second
        /// argument and the fact that the e-geodesic traces a straight line in natural
        /// parameter space (Amari & Nagaoka, 2000, Ch. 3).
        #[test]
        fn e_geodesic_kl_monotone(
            p in positive_simplex_vec(5),
            q in positive_simplex_vec(5),
        ) {
            let mut prev_kl = 0.0_f64; // KL(p || gamma(0)) = KL(p || p) = 0
            for step in 1..=10 {
                let t = step as f64 / 10.0;
                let gamma_t = e_geodesic(&p, &q, t, TOL).unwrap();
                let kl_t = logp::kl_divergence(&p, &gamma_t, TOL).unwrap();
                prop_assert!(kl_t >= prev_kl - 1e-10,
                    "KL not monotone: KL(p, gamma({t}))={kl_t} < prev={prev_kl}");
                prev_kl = kl_t;
            }
        }
    }

    proptest! {
        /// For uniform p = [1/n, ..., 1/n], natural_gradient(p, g) = g / n.
        #[test]
        fn natural_gradient_on_uniform_equals_euclidean_scaled(
            g in prop::collection::vec(-10.0f64..10.0, 5),
        ) {
            let n = g.len();
            let p: Vec<f64> = vec![1.0 / n as f64; n];
            let ng = natural_gradient(&p, &g).unwrap();

            for i in 0..n {
                let expected = g[i] / n as f64;
                prop_assert!((ng[i] - expected).abs() < 1e-14,
                    "natural_gradient[{i}]={} != g[{i}]/n={}", ng[i], expected);
            }
        }
    }

    #[test]
    fn natural_gradient_length_mismatch_errors() {
        let p = vec![0.5, 0.5];
        let g = vec![1.0, 2.0, 3.0];
        let result = natural_gradient(&p, &g);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("length mismatch"),
            "expected length mismatch error, got: {err_msg}"
        );
    }

    #[test]
    fn fisher_rao_geodesic_invalid_t_errors() {
        let p = vec![0.5, 0.3, 0.2];
        let q = vec![0.2, 0.3, 0.5];

        let neg = fisher_rao_geodesic(&p, &q, -0.1, 1e-12);
        assert!(neg.is_err());
        assert!(neg.unwrap_err().to_string().contains("outside [0, 1]"));

        let over = fisher_rao_geodesic(&p, &q, 1.5, 1e-12);
        assert!(over.is_err());
        assert!(over.unwrap_err().to_string().contains("outside [0, 1]"));
    }

    #[test]
    fn alpha_geodesic_subsumes_m_and_e() {
        let p = [0.5, 0.3, 0.2];
        let q = [0.2, 0.5, 0.3];
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let a = alpha_geodesic(&p, &q, t, -1.0, 1e-12).unwrap();
            let m = m_geodesic(&p, &q, t, 1e-12).unwrap();
            for (x, y) in a.iter().zip(m.iter()) {
                assert!((x - y).abs() < 1e-9, "alpha=-1 != m_geodesic at t={t}");
            }
            let a1 = alpha_geodesic(&p, &q, t, 1.0, 1e-12).unwrap();
            let e = e_geodesic(&p, &q, t, 1e-12).unwrap();
            for (x, y) in a1.iter().zip(e.iter()) {
                assert!((x - y).abs() < 1e-9, "alpha=+1 != e_geodesic at t={t}");
            }
        }
    }

    #[test]
    fn alpha_geodesic_endpoints_and_simplex() {
        let p = [0.5, 0.3, 0.2];
        let q = [0.2, 0.5, 0.3];
        for &alpha in &[-1.0, -0.5, 0.0, 0.5, 1.0] {
            let g0 = alpha_geodesic(&p, &q, 0.0, alpha, 1e-12).unwrap();
            let g1 = alpha_geodesic(&p, &q, 1.0, alpha, 1e-12).unwrap();
            for (x, y) in g0.iter().zip(p.iter()) {
                assert!((x - y).abs() < 1e-9);
            }
            for (x, y) in g1.iter().zip(q.iter()) {
                assert!((x - y).abs() < 1e-9);
            }
            for &t in &[0.3, 0.5] {
                let g = alpha_geodesic(&p, &q, t, alpha, 1e-12).unwrap();
                let s: f64 = g.iter().sum();
                assert!((s - 1.0).abs() < 1e-9, "alpha={alpha} t={t} sum={s}");
            }
        }
    }

    #[test]
    fn alpha_geodesic_negative_exponent_needs_positive() {
        let p = [0.5, 0.0, 0.5];
        let q = [0.2, 0.5, 0.3];
        assert!(alpha_geodesic(&p, &q, 0.5, 2.0, 1e-12).is_err());
        assert!(alpha_geodesic(&p, &q, 0.5, 0.0, 1e-12).is_ok());
    }
}
