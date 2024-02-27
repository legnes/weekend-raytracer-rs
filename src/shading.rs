use super::hit::HitRecord;
use super::ray::Ray;
use super::vec::Vec3;
use rand::Rng;
use std::f64::consts::PI;

#[allow(dead_code)]
pub enum DiffuseModel {
    Lambert,
    Lambert3,
    Uniform,
}

pub fn diffuse(hit: HitRecord, method: DiffuseModel) -> Ray {
    match method {
        // Diffuse lambertian --
        // Reflect in the direction of the normal modulated by a random vector on the unit sphere.
        // This scales by cos(theta) where theta is angle from the normal
        DiffuseModel::Lambert => {
            Ray::new(hit.position, hit.normal + random_on_unit_sphere())
        }
        // Approximate lambertian --
        // Reflect in the direction of the normal modulated by a random vector IN the unit sphere.
        // This scales by cos^3(theta), so it will be a bit darker
        // since more rays are scattered towards the normal
        // SE TODO: Not quite sure why this would mean darker...
        DiffuseModel::Lambert3 => {
            Ray::new(hit.position, hit.normal + random_in_unit_sphere())
        }
        // Uniform diffuse --
        // A uniform scatter direction for all angles away from the hit point,
        // with no dependence on the angle from the normal
        DiffuseModel::Uniform => {
            let mut dir = random_in_unit_sphere();
            if hit.normal.dot(dir) <= 0.0 {
                // In the opposite hemisphere as the normal
                dir *= -1.0;
            }
            Ray::new(hit.position, dir)
        }
    }
}

fn random_on_unit_sphere() -> Vec3 {
    let mut rng = rand::thread_rng();

    let phi = rng.gen_range(0.0..(2.0 * PI));
    // Polar angle calc from https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
    let theta = f64::acos(rng.gen_range(-1.0..1.0));

    Vec3::from_spherical(1.0, phi, theta)
}

fn random_in_unit_sphere() -> Vec3 {
    let mut rng = rand::thread_rng();

    let r = f64::powf(rng.gen(), 1.0 / 3.0);
    let phi = rng.gen_range(0.0..(2.0 * PI));
    let theta = f64::acos(rng.gen_range(-1.0..1.0));

    Vec3::from_spherical(r, phi, theta)
}
