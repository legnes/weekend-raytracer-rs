use super::ray::Ray;
use super::vec::{Color, Point3, Vec3};

pub trait Light: Send + Sync {
    fn color(&self) -> Color;

    // Ray from the surface position towards the light source, normalized by convention
    fn get_shadow_ray(&self, surface_position: Point3) -> Ray;
}

pub struct DirectionalLight {
    direction: Vec3,
    color: Color,
}

impl DirectionalLight {
    pub fn new(direction: Vec3, color: Color) -> Self {
        Self {
            direction: direction.normalized(),
            color,
        }
    }
}

impl Light for DirectionalLight {
    fn color(&self) -> Color {
        self.color
    }

    fn get_shadow_ray(&self, surface_position: Point3) -> Ray {
        Ray::new(surface_position, -self.direction)
    }
}
