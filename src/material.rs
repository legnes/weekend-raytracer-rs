use super::hit::HitRecord;
use super::ray::Ray;
use super::vec::{Color, Vec3};
use rand::Rng;

pub trait Scatter: Send + Sync {
    fn scatter(&self, incident: &Ray, hit: &HitRecord) -> Option<(Color, Ray)>;
    // We are assuming to_light and hit.normal are normalized
    fn shade(&self, from_view: &Ray, to_light: &Ray, hit: &HitRecord) -> Color;
    fn casts_shadow(&self) -> bool;
}

trait DiffuseModel: Send + Sync {
    fn get_scatter_direction(hit: &HitRecord) -> Vec3;
    fn get_albedo(&self) -> Color;
}

impl<T> Scatter for T
where
    T: DiffuseModel + Send + Sync,
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

    // SE TODO: This is lambertian, need to figure out the other two distributions
    fn shade(&self, _from_view: &Ray, to_light: &Ray, hit: &HitRecord) -> Color {
        self.get_albedo() * hit.normal.dot(to_light.direction()).clamp(0.0, 1.0)
    }

    fn casts_shadow(&self) -> bool {
        true
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
        hit.normal + Vec3::random_on_unit_sphere()
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
        let reflected = incident.direction().reflect(hit.normal);
        let scattered = Ray::new(
            hit.position,
            reflected + self.fuzz * Vec3::random_in_unit_sphere(),
        );

        // Make sure the fuzz has not put us inside the surface
        if scattered.direction().dot(hit.normal) > 0.0 {
            Some((self.albedo, scattered))
        } else {
            None
        }
    }

    fn shade(&self, from_view: &Ray, to_light: &Ray, hit: &HitRecord) -> Color {
        let view_reflected = from_view.direction().normalized().reflect(hit.normal);
        // Phong specular
        self.albedo
            * view_reflected
                .dot(to_light.direction())
                .clamp(0.0, 1.0)
                .powf(16.0 * (1.0 - self.fuzz))
    }

    fn casts_shadow(&self) -> bool {
        true
    }
}

pub struct Dielectric {
    refractive_index: f64,
}

impl Dielectric {
    #[allow(dead_code)]
    pub fn new(refractive_index: f64) -> Self {
        Self { refractive_index }
    }

    fn reflectance(cos: f64, refractive_index_1: f64, refractive_index_2: f64) -> f64 {
        // Schlick approximation
        let r0 = ((refractive_index_1 - refractive_index_2)
            / (refractive_index_1 + refractive_index_2))
            .powi(2);
        r0 + (1.0 - r0) * (1.0 - cos).powi(5)
    }
}

impl Scatter for Dielectric {
    fn scatter(&self, incident: &Ray, hit: &HitRecord) -> Option<(Color, Ray)> {
        // SE TODO: Looks like we are assuming this is always interacting with air?
        let other_refractive_index = 1.0;
        let refractive_ratio = if hit.front_face {
            other_refractive_index / self.refractive_index
        } else {
            self.refractive_index / other_refractive_index
        };

        let unit_incident = incident.direction().normalized();
        // SE TODO: Why do we need min() if already normalized?
        let cos_theta = (-unit_incident).dot(hit.normal).min(1.0);
        let sin_theta = (1.0 - cos_theta.powi(2)).sqrt();

        // Total internal reflection
        // If going from low to high index of refraction, there may not be a solution to Snell's law
        // when the incident angle is low enough (sin(theta) is large, so that ratio * sin(theta) > 1,
        // since sin(theta') can't be > 1)
        let cannot_refract = refractive_ratio * sin_theta > 1.0;

        // Fresnel reflection
        let mut rng = rand::thread_rng();
        // SE TODO: The book uses the ratio and 1.0 here? Check math and make sure it doesn't matter
        let will_reflect = rng.gen::<f64>()
            < Self::reflectance(cos_theta, self.refractive_index, other_refractive_index);

        let out_direction = if cannot_refract || will_reflect {
            // Reflect
            unit_incident.reflect(hit.normal)
        } else {
            // Refract
            unit_incident.refract(hit.normal, refractive_ratio)
        };

        let scattered = Ray::new(hit.position, out_direction);

        Some((Color::one(), scattered))
    }

    fn shade(&self, from_view: &Ray, to_light: &Ray, hit: &HitRecord) -> Color {
        let view_reflected = from_view.direction().normalized().reflect(hit.normal);
        Color::one()
            * view_reflected
                .dot(to_light.direction())
                .clamp(0.0, 1.0)
                .powi(16)
    }

    fn casts_shadow(&self) -> bool {
        // SE TODO: Need to make this conditional on stuff, factor in refraction, etc.
        false
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
        hit.normal + Vec3::random_in_unit_sphere()
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
        let dir = Vec3::random_in_unit_sphere();
        if hit.normal.dot(dir) <= 0.0 {
            // In the opposite hemisphere as the normal
            -dir
        } else {
            dir
        }
    }
}
