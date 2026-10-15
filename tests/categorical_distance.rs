use infogeom::rao_distance_categorical as distance;

fn bernoulli(x: f64) -> [f64; 2] {
    [x, 1.0 - x]
}

// The Bernoulli arc coordinate is 2 asin(sqrt(x)), independent of BC/chords.
fn angular_distance(x: f64, y: f64) -> f64 {
    2.0 * (x.sqrt().asin() - y.sqrt().asin()).abs()
}

#[test]
fn nearby_triangle_preserves_bernoulli_arc_lengths() {
    for (tol, step) in [(1e-3, 0.08), (1e-9, 0.0001), (1e-12, 0.000003)] {
        let (p, q, r) = (
            bernoulli(0.5),
            bernoulli(0.5 + step),
            bernoulli(0.5 + 2.0 * step),
        );
        let pq = distance(&p, &q, tol).unwrap();
        let qr = distance(&q, &r, tol).unwrap();
        let pr = distance(&p, &r, tol).unwrap();
        assert!(pq > 0.0 && qr > 0.0);
        assert!(pr <= pq + qr + 2e-15);
        for (x, y, actual) in [(p[0], q[0], pq), (q[0], r[0], qr), (p[0], r[0], pr)] {
            assert!((actual - angular_distance(x, y)).abs() < 2e-15);
        }
    }
}

#[test]
fn tiny_boundary_and_adjacent_probabilities_remain_visible() {
    let tiny = f64::from_bits(1);
    let actual = distance(&[1.0, 0.0], &[1.0, tiny], 1e-12).unwrap();
    let expected = 2.0 * tiny.sqrt();
    assert!(actual > 0.0);
    assert!((actual / expected - 1.0).abs() < 1e-14);
    let x = 0.5_f64;
    let y = f64::from_bits(x.to_bits() + 1);
    let actual = distance(&bernoulli(x), &bernoulli(y), 1e-12).unwrap();
    // At one half, the Bernoulli arc derivative is exactly 2.
    assert!((actual / (2.0 * (y - x)) - 1.0).abs() < 1e-14);
}

#[test]
fn normalization_preserves_root_representable_boundary_mass() {
    let tiny = f64::from_bits(1);
    for mass in [2.0_f64, 1e308] {
        let actual = distance(&[tiny, mass], &[0.0, mass], mass).unwrap();
        let expected = 2.0 * (tiny.sqrt() / mass.sqrt());
        assert!(actual > 0.0);
        assert!((actual - expected).abs() <= 4.0 * f64::from_bits(1));
    }
    // Different totals exercise the normalization branch with an underflowing ratio.
    let actual = distance(&[tiny, 2.0], &[0.0, 1.0], 1.0).unwrap();
    assert!((actual / (2.0 * tiny.sqrt() / 2.0_f64.sqrt()) - 1.0).abs() < 1e-14);
    let actual = distance(&[3.0 * tiny, 2.0], &[tiny, 1.0], 1.0).unwrap();
    let expected = 2.0 * ((3.0 * tiny).sqrt() / 2.0_f64.sqrt() - tiny.sqrt());
    assert!((actual / expected - 1.0).abs() < 1e-14);
}

#[test]
fn validation_tolerance_does_not_change_admitted_distances() {
    let p = [0.5, 0.5];
    let q = [0.500003, 0.499997];
    let expected = distance(&p, &q, 0.0).unwrap();
    for tol in [1e-12, 1e-9, 1e-3, 1.0] {
        assert_eq!(distance(&p, &q, tol).unwrap(), expected);
    }
}

#[test]
fn accepted_positive_mass_is_normalized() {
    let p = [0.25, 0.75];
    let q = [0.75, 0.25];
    assert_eq!(distance(&p, &[0.5, 1.5], 1.0).unwrap(), 0.0);
    let expected = distance(&p, &q, 1e-12).unwrap();
    let scaled = [q[0] * (1.0 + 1e-8), q[1] * (1.0 + 1e-8)];
    assert!((distance(&p, &scaled, 2e-8).unwrap() - expected).abs() < 1e-14);
    assert!(distance(&p, &scaled, 1e-12).is_err());
    assert_eq!(distance(&[0.99], &[1.01], 0.02).unwrap(), 0.0);
}

#[test]
fn disjoint_support_attains_the_exact_range_endpoint() {
    for (p, q) in [
        (vec![1.0, 0.0], vec![0.0, 1.0]),
        (vec![0.25, 0.75, 0.0], vec![0.0, 0.0, 1.0]),
    ] {
        let actual = distance(&p, &q, 1e-12).unwrap();
        assert!(
            actual <= std::f64::consts::PI,
            "distance={actual:.17e}, pi={:.17e}",
            std::f64::consts::PI
        );
        assert_eq!(actual, std::f64::consts::PI);
    }
}

#[test]
fn self_symmetry_boundary_and_permutation() {
    let p = [0.25, 0.0, 0.75];
    let q = [0.75, 0.0, 0.25];
    assert_eq!(distance(&p, &p, 1e-12).unwrap(), 0.0);
    let d = distance(&p, &q, 1e-12).unwrap();
    assert_eq!(distance(&q, &p, 1e-12).unwrap(), d);
    assert_eq!(
        distance(&[0.0, 0.75, 0.25], &[0.0, 0.25, 0.75], 1e-12).unwrap(),
        d
    );
    assert!((d - std::f64::consts::PI / 3.0).abs() < 1e-14);
    assert!(
        (distance(&[1.0, 0.0], &[0.0, 1.0], 1.0).unwrap() - std::f64::consts::PI).abs() < 1e-14
    );
}

#[test]
fn invalid_tolerance_mass_and_existing_input_errors() {
    for tol in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
        assert!(matches!(
            distance(&[1.0], &[1.0], tol),
            Err(infogeom::Error::Logp(logp::Error::Domain(_)))
        ));
    }
    assert!(matches!(
        distance(&[0.0], &[1.0], 1.0),
        Err(infogeom::Error::Logp(logp::Error::Domain(_)))
    ));
    assert!(matches!(
        distance(&[], &[1.0], 1e-12),
        Err(infogeom::Error::Logp(logp::Error::LengthMismatch(0, 1)))
    ));
    assert!(matches!(
        distance(&[], &[], 1e-12),
        Err(infogeom::Error::Logp(logp::Error::Empty))
    ));
    assert!(matches!(
        distance(&[f64::NAN], &[1.0], 1e-12),
        Err(infogeom::Error::Logp(logp::Error::NonFinite { .. }))
    ));
    assert!(matches!(
        distance(&[-1.0, 2.0], &[0.0, 1.0], 1e-12),
        Err(infogeom::Error::Logp(logp::Error::Negative { .. }))
    ));
    assert!(matches!(
        distance(&[0.5], &[1.0], 1e-12),
        Err(infogeom::Error::Logp(logp::Error::NotNormalized { .. }))
    ));
    // Keep the existing simplex sum error when a finite-entry total overflows.
    assert!(matches!(
        distance(&[f64::MAX, f64::MAX], &[0.5, 0.5], f64::MAX),
        Err(infogeom::Error::Logp(logp::Error::NotNormalized { .. }))
    ));
}
