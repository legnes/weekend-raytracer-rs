use super::bvh::Bvh;
use super::material::Scatter;
use super::ray::Ray;
use super::vec::{Point3, Vec3};
use std::sync::Arc;

pub struct HitRecord {
    pub t: f64,
    pub position: Point3,
    pub normal: Vec3,
    pub front_face: bool,
    pub material: Arc<dyn Scatter>,
}

impl HitRecord {
    pub fn set_face_normal(&mut self, ray: &Ray, outward_normal: Vec3) -> () {
        // Point normal against ray for faster shading calc (skip dot prod)
        // As a result we have to track if this is front/back face of surface
        // Could instead store normal as always outwards and use dot prod during shading
        self.front_face = ray.direction().dot(outward_normal) < 0.0;
        self.normal = if self.front_face {
            outward_normal
        } else {
            -outward_normal
        }
    }
}

pub trait Hit : Send + Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

pub trait Bounded : Send + Sync {
    fn centroid(&self) -> Point3;
    fn aabb_min(&self) -> Point3;
    fn aabb_max(&self) -> Point3;
}

pub trait Primitive : Hit + Bounded {}

pub type World = Vec<Box<dyn Primitive>>;

pub struct Scene {
    primitives: World,
    bvh: Bvh,
}

impl Scene {
    pub fn new(primitives: World) -> Self {
        let bvh = Bvh::new(&primitives);
        Self { primitives, bvh }
    }

    #[allow(dead_code)]
    fn hit_slow(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest_t = t_max;
        let mut closest_record = None;

        for primitive in &self.primitives {
            if let Some(record) = primitive.hit(ray, t_min, closest_t) {
                closest_t = record.t;
                closest_record = Some(record);
            }
        }

        closest_record
    }
}

impl Hit for Scene {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        // self.hit_slow(ray, t_min, t_max)
        self.bvh.hit(ray, t_min, t_max, &self.primitives)
    }
}
