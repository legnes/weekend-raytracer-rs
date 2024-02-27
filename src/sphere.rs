use super::hit::{Hit, HitRecord};
use super::material::Scatter;
use super::ray::Ray;
use super::vec::{Point3, Vec3};
use std::rc::Rc;

pub struct Sphere {
    center: Point3,
    radius: f64,
    material: Rc<dyn Scatter>,
}

impl Sphere {
    pub fn new(center: Point3, radius: f64, material: Rc<dyn Scatter>) -> Self {
        Self {
            center,
            radius,
            material,
        }
    }
}

impl Hit for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin() - self.center;
        // Quadratic formula (simplified)
        let a = ray.direction().length().powi(2);
        let half_b = oc.dot(ray.direction());
        let c = oc.length().powi(2) - self.radius.powi(2);

        let discriminant = half_b.powi(2) - a * c;
        if discriminant < 0.0 {
            return None;
        }

        // Find the nearest root that lies in the acceptable range
        let sqrt_d = discriminant.sqrt();
        // This is the smallest root from the (simplified) quadratic formula
        let mut root = (-half_b - sqrt_d) / a;
        if root < t_min || root > t_max {
            root = (-half_b + sqrt_d) / a;
            if root < t_min || root > t_max {
                return None;
            }
        }

        let position = ray.at(root);
        let mut hit = HitRecord {
            t: root,
            position,
            normal: Vec3::zero(),
            front_face: false,
            material: self.material.clone(),
        };

        let outward_normal = (position - self.center) / self.radius;
        hit.set_face_normal(ray, outward_normal);

        Some(hit)
    }
}
