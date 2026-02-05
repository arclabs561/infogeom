use infogeom::{hellinger, rao_distance_categorical};

fn main() {
    let p = [0.70, 0.20, 0.10];
    let q = [0.10, 0.20, 0.70];

    let d_rao = rao_distance_categorical(&p, &q, 1e-12).unwrap();
    let d_hel = hellinger(&p, &q, 1e-12).unwrap();

    assert!(d_rao >= 0.0);
    assert!((0.0..=1.0).contains(&d_hel));

    println!("Rao distance (radians): {:.6}", d_rao);
    println!("Hellinger distance:     {:.6}", d_hel);
}

