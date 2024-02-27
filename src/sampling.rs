use super::vec::Vec3;
use rand::Rng;
use std::f64::consts::PI;

pub fn random_on_unit_sphere() -> Vec3 {
    let mut rng = rand::thread_rng();

    let phi = rng.gen_range(0.0..(2.0 * PI));
    // Polar angle calc from https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
    let theta = f64::acos(rng.gen_range(-1.0..1.0));

    Vec3::from_spherical(1.0, phi, theta)
}

pub fn random_in_unit_sphere() -> Vec3 {
    let mut rng = rand::thread_rng();

    let r = f64::powf(rng.gen(), 1.0 / 3.0);
    let phi = rng.gen_range(0.0..(2.0 * PI));
    let theta = f64::acos(rng.gen_range(-1.0..1.0));

    Vec3::from_spherical(r, phi, theta)
}
