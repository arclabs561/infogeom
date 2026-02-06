//! `infogeom`: information geometry primitives.
//!
//! This crate provides small, policy-free building blocks for geometry on probability
//! distributions.
//!
//! Today it focuses on the probability simplex (categorical distributions):
//! - Fisher–Rao / Rao distance
//! - Hellinger distance
//!
//! `logp` provides divergence/entropy functionals and simplex validation; `infogeom` builds
//! geometry on top.
//!
//! Related:
//! - `unseen` estimates properties from samples (unseen regime) when you don’t have an explicit
//!   distribution yet.
//!
//! Reference orientation:
//! - Frank Nielsen’s divergence/IG portal: https://franknielsen.github.io/IG/index.html
//!
//! ## Quick example
//!
//! ```rust
//! use infogeom::{rao_distance_categorical, hellinger};
//!
//! let p = [0.7, 0.2, 0.1];
//! let q = [0.1, 0.2, 0.7];
//!
//! let rao = rao_distance_categorical(&p, &q, 1e-12).unwrap();
//! let hel = hellinger(&p, &q, 1e-12).unwrap();
//!
//! assert!(rao >= 0.0);
//! assert!((0.0..=1.0).contains(&hel));
//! ```

#![forbid(unsafe_code)]

use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    Logp(#[from] logp::Error),

    #[error("domain error: {0}")]
    Domain(&'static str),
}

pub type Result<T> = core::result::Result<T, Error>;

/// Fisher–Rao / Rao distance between categorical distributions `p` and `q` (in radians).
///
/// For categorical distributions, the Fisher–Rao manifold can be isometrically embedded on a
/// sphere via the map \(p \mapsto 2\sqrt{p}\). Under this embedding, Rao distance is twice the
/// spherical angle:
///
/// \( d_{FR}(p,q) = 2 \arccos\left(\sum_i \sqrt{p_i q_i}\right) \).
///
/// Notes:
/// - The inner term is the Bhattacharyya coefficient \(BC(p,q)\in[0,1]\).
/// - We validate both inputs are on the simplex within `tol`.
pub fn rao_distance_categorical(p: &[f64], q: &[f64], tol: f64) -> Result<f64> {
    // Delegate simplex checks to `logp` and reuse its Bhattacharyya coefficient.
    let mut bc = logp::bhattacharyya_coeff(p, q, tol)?;

    // Numerical handling:
    // - rounding can push slightly outside [0,1]
    // - when p==q, bc should be 1 but may be 1-ε which makes acos nonzero
    bc = bc.clamp(0.0, 1.0);
    if (1.0 - bc).abs() <= 10.0 * tol {
        bc = 1.0;
    }
    Ok(2.0 * bc.acos())
}

/// Hellinger distance \(H(p,q)\) induced by the simplex sphere embedding.
///
/// For categorical distributions, \(H^2(p,q) = 1 - BC(p,q)\), where \(BC\) is the Bhattacharyya
/// coefficient. This is a bounded metric on the simplex with range \([0, 1]\).
pub fn hellinger(p: &[f64], q: &[f64], tol: f64) -> Result<f64> {
    Ok(logp::hellinger(p, q, tol)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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

    proptest! {
        #[test]
        fn rao_is_symmetric(p in simplex_vec(8), q in simplex_vec(8)) {
            let d1 = rao_distance_categorical(&p, &q, 1e-9).unwrap();
            let d2 = rao_distance_categorical(&q, &p, 1e-9).unwrap();
            prop_assert!((d1 - d2).abs() < 1e-12);
        }

        #[test]
        fn rao_self_is_zero(p in simplex_vec(12)) {
            let d = rao_distance_categorical(&p, &p, 1e-9).unwrap();
            prop_assert!(d.abs() < 1e-12);
        }

        #[test]
        fn rao_is_bounded(p in simplex_vec(10), q in simplex_vec(10)) {
            let d = rao_distance_categorical(&p, &q, 1e-9).unwrap();
            // Since bc ∈ [0,1], acos(bc) ∈ [0, π/2], so 2*acos(bc) ∈ [0, π].
            prop_assert!(d >= -1e-12);
            prop_assert!(d <= core::f64::consts::PI + 1e-12);
        }

        #[test]
        fn hellinger_is_symmetric(p in simplex_vec(8), q in simplex_vec(8)) {
            let h1 = hellinger(&p, &q, 1e-9).unwrap();
            let h2 = hellinger(&q, &p, 1e-9).unwrap();
            prop_assert!((h1 - h2).abs() < 1e-12);
        }

        #[test]
        fn hellinger_is_bounded(p in simplex_vec(8), q in simplex_vec(8)) {
            let h = hellinger(&p, &q, 1e-9).unwrap();
            prop_assert!(h >= -1e-12);
            prop_assert!(h <= 1.0 + 1e-12);
        }

        #[test]
        fn rao_matches_bhattacharyya_formula(p in simplex_vec(10), q in simplex_vec(10)) {
            let tol = 1e-9;
            let rao = rao_distance_categorical(&p, &q, tol).unwrap();
            let mut bc = logp::bhattacharyya_coeff(&p, &q, tol).unwrap();
            bc = bc.clamp(0.0, 1.0);
            if (1.0 - bc).abs() <= 10.0 * tol {
                bc = 1.0;
            }
            let expected = 2.0 * bc.acos();
            prop_assert!((rao - expected).abs() < 1e-12, "rao={rao} expected={expected} bc={bc}");
        }
    }
}

