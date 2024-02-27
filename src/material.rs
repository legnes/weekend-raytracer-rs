use crate::sampling::random_in_unit_sphere;

use super::hit::HitRecord;
use super::ray::Ray;
use super::sampling;
use super::vec::{Color, Vec3};

pub trait Scatter {
    fn scatter(&self, incident: &Ray, hit: &HitRecord) -> Option<(Color, Ray)>;
}

trait DiffuseModel {
    fn get_scatter_direction(hit: &HitRecord) -> Vec3;
    fn get_albedo(&self) -> Color;
}

impl<T> Scatter for T
where
    T: DiffuseModel,
{
    // Note we could just as well only scatter with some probability p and have attenuation be albedo/p
    fn scatter(&self, _incident: &Ray, hit: &HitRecord) -> Option<(Color, Ray)> {
        let mut scatter_direction = Self::get_scatter_direction(hit);

        if scatter_direction.near_zero() {
            // Degenerate case -- we have sampled the direction exactly opposite of the normal,
            // causing a ray with little/no magnitude. Just scatter along normal.
            scatter_direction = hit.normal;
        }

        Some((self.get_albedo(), Ray::new(hit.position, scatter_direction)))
    }
}

pub struct Lambertian {
    albedo: Color,
}

impl Lambertian {
    #[allow(dead_code)]
    pub fn new(albedo: Color) -> Self {
        Self { albedo }
    }
}

impl DiffuseModel for Lambertian {
    fn get_albedo(&self) -> Color {
        self.albedo
    }

    fn get_scatter_direction(hit: &HitRecord) -> Vec3 {
        // Diffuse lambertian --
        // Reflect in the direction of the normal modulated by a random vector on the unit sphere.
        // This scales by cos(theta) where theta is angle from the normal
        hit.normal + sampling::random_on_unit_sphere()
    }
}

pub struct Metal {
    albedo: Color,
    fuzz: f64,
}

impl Metal {
    #[allow(dead_code)]
    pub fn new(albedo: Color, fuzz: f64) -> Self {
        Self { albedo, fuzz }
    }
}

impl Scatter for Metal {
    fn scatter(&self, incident: &Ray, hit: &HitRecord) -> Option<(Color, Ray)> {
        let scatter_direction =
            incident.direction().reflect(hit.normal) + self.fuzz * random_in_unit_sphere();
        let reflection = Ray::new(hit.position, scatter_direction);

        // Make sure the fuzz has not put us inside the surface
        if reflection.direction().dot(hit.normal) > 0.0 {
            Some((self.albedo, reflection))
        } else {
            None
        }
    }
}

//////////////////////////////////////
//
// Other models
//
//////////////////////////////////////
pub struct ApproximateLambertian {
    albedo: Color,
}

impl ApproximateLambertian {
    #[allow(dead_code)]
    pub fn new(albedo: Color) -> Self {
        Self { albedo }
    }
}

impl DiffuseModel for ApproximateLambertian {
    fn get_albedo(&self) -> Color {
        self.albedo
    }

    fn get_scatter_direction(hit: &HitRecord) -> Vec3 {
        // Approximate Lambertian --
        // Reflect in the direction of the normal modulated by a random vector IN the unit sphere.
        // This scales by cos^3(theta), so it will be a bit darker
        // since more rays are scattered towards the normal
        // SE TODO: Not quite sure why this would mean darker...
        hit.normal + sampling::random_in_unit_sphere()
    }
}

pub struct UniformDiffuse {
    albedo: Color,
}

impl UniformDiffuse {
    #[allow(dead_code)]
    pub fn new(albedo: Color) -> Self {
        Self { albedo }
    }
}

impl DiffuseModel for UniformDiffuse {
    fn get_albedo(&self) -> Color {
        self.albedo
    }

    fn get_scatter_direction(hit: &HitRecord) -> Vec3 {
        // Uniform Diffuse --
        // A uniform scatter direction for all angles away from the hit point,
        // with no dependence on the angle from the normal
        let dir = sampling::random_in_unit_sphere();
        if hit.normal.dot(dir) <= 0.0 {
            // In the opposite hemisphere as the normal
            (-1.0) * dir
        } else {
            dir
        }
    }
}
