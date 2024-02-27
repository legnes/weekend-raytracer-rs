use super::vec::{Point3, Vec3};

pub struct Ray {
    origin: Point3,
    direction: Vec3,
}

impl Ray {
    pub fn new(origin: Point3, direction: Vec3) -> Self {
        Self { origin, direction }
    }

    // SE TODO: Why &self here but self for Vec3.x() etc.?
    // Answer: because vec3 implements Copy, so is stored on the stack
    // and rust always passes it by copy
    // SE TODO: Why getter functions for these? Rust enforces immutability...?
    // SE TODO: Is it idiomatic to name the member and the getter differently?
    pub fn origin(&self) -> Point3 {
        self.origin
    }

    pub fn direction(&self) -> Vec3 {
        self.direction
    }

    pub fn at(&self, t: f64) -> Point3 {
        self.origin + t * self.direction
    }
}
