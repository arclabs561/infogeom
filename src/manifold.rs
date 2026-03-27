//! [`skel::Manifold`] implementation for the Fisher-Rao simplex.
//!
//! The categorical simplex under the Fisher-Rao metric embeds isometrically
//! into the positive orthant of the unit sphere via `theta_i = sqrt(p_i)`.
//! All Riemannian operations (exp, log, transport) reduce to standard sphere
//! operations in this embedding, then pull back to the simplex by squaring.
//!
//! Requires the `manifold` feature.

use ndarray::{Array1, ArrayView1};
use skel::Manifold;

/// Fisher-Rao Riemannian structure on the probability simplex.
///
/// Points are probability vectors `p = [p_1, ..., p_n]` with `p_i >= 0`,
/// `sum(p_i) = 1`. Tangent vectors at `p` satisfy `sum(v_i / (2*sqrt(p_i))) = 0`
/// in the sphere embedding.
///
/// The tolerance parameter controls the threshold below which two points are
/// considered identical (avoiding division by zero in log_map).
#[derive(Debug, Clone)]
pub struct FisherRaoSimplex {
    tol: f64,
}

impl FisherRaoSimplex {
    /// Create a new Fisher-Rao simplex manifold with the given tolerance.
    pub fn new(tol: f64) -> Self {
        Self { tol }
    }
}

impl Default for FisherRaoSimplex {
    fn default() -> Self {
        Self { tol: 1e-12 }
    }
}

// Helpers: simplex <-> sphere embedding.
//
// Embedding:  theta_i = sqrt(p_i)         (simplex -> sphere)
// Pullback:   p_i = theta_i^2             (sphere -> simplex)
//
// Tangent lift: if v is a tangent vector at p on the simplex, the
// corresponding tangent on the sphere is:
//   u_i = v_i / (2 * sqrt(p_i))
//
// Tangent drop: the inverse:
//   v_i = 2 * sqrt(p_i) * u_i

fn to_sphere(p: &ArrayView1<f64>) -> Array1<f64> {
    p.mapv(|pi| pi.max(0.0).sqrt())
}

fn to_simplex(theta: &Array1<f64>) -> Array1<f64> {
    let sq = theta.mapv(|t| t * t);
    let s: f64 = sq.sum();
    if s > 0.0 {
        sq / s
    } else {
        sq
    }
}

fn tangent_lift(p: &ArrayView1<f64>, v: &ArrayView1<f64>) -> Array1<f64> {
    Array1::from_iter(p.iter().zip(v.iter()).map(|(&pi, &vi)| {
        let sq = pi.max(1e-300).sqrt();
        vi / (2.0 * sq)
    }))
}

fn tangent_drop(p: &ArrayView1<f64>, u: &Array1<f64>) -> Array1<f64> {
    Array1::from_iter(
        p.iter()
            .zip(u.iter())
            .map(|(&pi, &ui)| 2.0 * pi.max(0.0).sqrt() * ui),
    )
}

impl Manifold for FisherRaoSimplex {
    fn exp_map(&self, x: &ArrayView1<f64>, v: &ArrayView1<f64>) -> Array1<f64> {
        let theta = to_sphere(x);
        let u = tangent_lift(x, v);

        // Geodesic on the sphere: slerp from theta in direction u.
        let norm = u.dot(&u).sqrt();
        if norm < self.tol {
            return x.to_owned();
        }

        let result = &theta * norm.cos() + &u * (norm.sin() / norm);
        to_simplex(&result)
    }

    fn log_map(&self, x: &ArrayView1<f64>, y: &ArrayView1<f64>) -> Array1<f64> {
        let theta_x = to_sphere(x);
        let theta_y = to_sphere(y);

        // Angle between the two points on the sphere.
        let dot = theta_x.dot(&theta_y).clamp(-1.0, 1.0);
        let omega = dot.acos();

        if omega < self.tol {
            return Array1::zeros(x.len());
        }

        // Direction on the sphere: project theta_y onto tangent space at theta_x.
        let u_unnorm = &theta_y - &theta_x * dot;
        let u_norm = u_unnorm.dot(&u_unnorm).sqrt();
        if u_norm < self.tol {
            return Array1::zeros(x.len());
        }
        let u = &u_unnorm * (omega / u_norm);

        // Drop back to simplex tangent space.
        tangent_drop(x, &u)
    }

    fn parallel_transport(
        &self,
        x: &ArrayView1<f64>,
        y: &ArrayView1<f64>,
        v: &ArrayView1<f64>,
    ) -> Array1<f64> {
        let theta_x = to_sphere(x);
        let theta_y = to_sphere(y);
        let u = tangent_lift(x, v);

        // Sphere parallel transport: Schild's ladder / closed-form.
        // For unit sphere: PT_{x->y}(u) = u - (dot(u, log) / omega^2) * (log_x + log_y)
        // where log = log_x(y) on the sphere.
        let dot_xy = theta_x.dot(&theta_y).clamp(-1.0, 1.0);
        let omega = dot_xy.acos();

        if omega < self.tol {
            return v.to_owned();
        }

        // log_x(y) on the sphere (unnormalized direction, scaled by omega).
        let proj = &theta_y - &theta_x * dot_xy;
        let proj_norm = proj.dot(&proj).sqrt();
        if proj_norm < self.tol {
            return v.to_owned();
        }
        let e = &proj / proj_norm; // unit tangent direction

        // Sphere parallel transport formula:
        // PT(u) = u - (theta_x * sin(omega) + e * (1 - cos(omega))) * dot(e, u)
        //       - (e * sin(omega) - theta_x * (1 - cos(omega))) * dot(theta_x, u)
        // Simplified for unit vectors:
        let eu = e.dot(&u);
        let xu = theta_x.dot(&u);
        let transported = &u
            + &e * (-omega.sin() * xu + (omega.cos() - 1.0) * eu)
            + &theta_x * (-omega.sin() * eu - (omega.cos() - 1.0) * xu);

        tangent_drop(y, &transported)
    }

    fn project(&self, x: &ArrayView1<f64>) -> Array1<f64> {
        // Project onto the simplex: clamp negatives, renormalize.
        let clamped = x.mapv(|v| v.max(0.0));
        let s: f64 = clamped.sum();
        if s > 0.0 {
            clamped / s
        } else {
            // Fallback: uniform distribution.
            let n = x.len() as f64;
            Array1::from_elem(x.len(), 1.0 / n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use proptest::prelude::*;

    fn manifold() -> FisherRaoSimplex {
        FisherRaoSimplex::new(1e-12)
    }

    /// Generate a strictly positive probability vector as Array1.
    fn positive_simplex(len: usize) -> impl Strategy<Value = Array1<f64>> {
        prop::collection::vec(0.01f64..10.0, len).prop_map(|v| {
            let s: f64 = v.iter().sum();
            Array1::from_iter(v.iter().map(|&x| x / s))
        })
    }

    proptest! {
        /// Round-trip: exp(log(y)) = y.
        #[test]
        fn round_trip_exp_log(
            p in positive_simplex(5),
            q in positive_simplex(5),
        ) {
            let m = manifold();
            let v = m.log_map(&p.view(), &q.view());
            let q_recovered = m.exp_map(&p.view(), &v.view());

            for i in 0..5 {
                prop_assert!((q_recovered[i] - q[i]).abs() < 1e-8,
                    "round-trip failed at {i}: recovered={} expected={}", q_recovered[i], q[i]);
            }
        }

        /// Round-trip: log(exp(v)) = v.
        #[test]
        fn round_trip_log_exp(
            p in positive_simplex(5),
            q in positive_simplex(5),
        ) {
            let m = manifold();
            let v = m.log_map(&p.view(), &q.view());
            let q2 = m.exp_map(&p.view(), &v.view());
            let v2 = m.log_map(&p.view(), &q2.view());

            for i in 0..5 {
                prop_assert!((v2[i] - v[i]).abs() < 1e-8,
                    "log-exp round-trip failed at {i}: v2={} v={}", v2[i], v[i]);
            }
        }

        /// log_map norm equals Fisher-Rao distance.
        #[test]
        fn log_norm_equals_rao_distance(
            p in positive_simplex(6),
            q in positive_simplex(6),
        ) {
            let m = manifold();
            let v = m.log_map(&p.view(), &q.view());

            // Compute norm in Fisher metric: sum(v_i^2 / (4 * p_i))
            // This equals the sphere-side norm ||u||^2 = omega^2.
            let sphere_u = tangent_lift(&p.view(), &v.view());
            let log_norm = sphere_u.dot(&sphere_u).sqrt();

            let rao = crate::rao_distance_categorical(
                p.as_slice().unwrap(), q.as_slice().unwrap(), 1e-12
            ).unwrap();

            // rao_distance = 2*arccos(BC), log_norm = arccos(dot) on the sphere.
            // The sphere embedding halves the distance: d_sphere = d_FR / 2.
            // But our tangent lift maps to the full sphere distance.
            prop_assert!((log_norm - rao / 2.0).abs() < 1e-8,
                "norm mismatch: log_norm={log_norm} rao/2={}", rao / 2.0);
        }

        /// Parallel transport preserves norm.
        #[test]
        fn transport_preserves_norm(
            p in positive_simplex(5),
            q in positive_simplex(5),
            r in positive_simplex(5),
        ) {
            let m = manifold();
            let v = m.log_map(&p.view(), &r.view());
            let tv = m.parallel_transport(&p.view(), &q.view(), &v.view());

            // Compute Fisher norms at p and q.
            let u_p = tangent_lift(&p.view(), &v.view());
            let u_q = tangent_lift(&q.view(), &tv.view());
            let norm_p = u_p.dot(&u_p).sqrt();
            let norm_q = u_q.dot(&u_q).sqrt();

            prop_assert!((norm_p - norm_q).abs() < 1e-6,
                "transport changed norm: {norm_p} -> {norm_q}");
        }

        /// Self-transport is identity.
        #[test]
        fn self_transport_identity(
            p in positive_simplex(5),
            q in positive_simplex(5),
        ) {
            let m = manifold();
            let v = m.log_map(&p.view(), &q.view());
            let tv = m.parallel_transport(&p.view(), &p.view(), &v.view());

            for i in 0..5 {
                prop_assert!((tv[i] - v[i]).abs() < 1e-10,
                    "self-transport not identity at {i}: tv={} v={}", tv[i], v[i]);
            }
        }
    }

    #[test]
    fn exp_zero_is_identity() {
        let m = manifold();
        let p = array![0.3, 0.5, 0.2];
        let zero = Array1::zeros(3);
        let result = m.exp_map(&p.view(), &zero.view());
        for i in 0..3 {
            assert!((result[i] - p[i]).abs() < 1e-14);
        }
    }

    #[test]
    fn log_self_is_zero() {
        let m = manifold();
        let p = array![0.3, 0.5, 0.2];
        let v = m.log_map(&p.view(), &p.view());
        for i in 0..3 {
            assert!(v[i].abs() < 1e-14, "log(p,p)[{i}] = {}", v[i]);
        }
    }

    #[test]
    fn project_renormalizes() {
        let m = manifold();
        let x = array![0.3, 0.5, -0.1, 0.4];
        let p = m.project(&x.view());
        assert!(p.iter().all(|&v| v >= 0.0));
        assert!(((p.sum()) - 1.0).abs() < 1e-14);
    }
}
