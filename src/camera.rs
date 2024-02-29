use super::ray::Ray;
use super::vec::{Point3, Vec3};
use std::f64::consts::PI;

pub struct Camera {
    origin: Point3,
    focus_plane_horizontal: Vec3,
    focus_plane_vertical: Vec3,
    focus_plane_lower_left_corner: Point3,
    camera_x: Vec3,
    camera_y: Vec3,
    lens_radius: f64,
}

impl Camera {
    pub fn new(
        origin: Point3,
        look_target: Point3,
        world_up: Vec3,
        vertical_fov_degrees: f64,
        aspect_ratio: f64,
        aperture: f64,
        focus_distance: f64,
    ) -> Self {
        // Viewport:
        // Doesn't really represent anything in the camera
        const FOCAL_LENGTH: f64 = 1.0;
        let fov = vertical_fov_degrees * (PI / 180.0);
        let viewport_height = 2.0 * FOCAL_LENGTH * (fov / 2.0).tan();
        let viewport_width = aspect_ratio * viewport_height;

        // Local coordinate system:
        // Camera faces along -z
        let camera_z = (origin - look_target).normalized();
        let camera_x = world_up.cross(camera_z).normalized();
        let camera_y = camera_z.cross(camera_x);

        // Focus plane:
        // Virtual film plane in the world, where rays from anywhere on the lens pass through
        // the same world point for a given uv, and all objects are in focus
        // SE TODO: I could be wrong about this math but I think it is using similar triangles
        // to get the dimensions of the focus plane from the dims of the viewport
        let focus_plane_horizontal = (focus_distance / FOCAL_LENGTH) * viewport_width * camera_x;
        let focus_plane_vertical = (focus_distance / FOCAL_LENGTH) * viewport_height * camera_y;
        let focus_plane_lower_left_corner = origin
            - focus_plane_horizontal / 2.0
            - focus_plane_vertical / 2.0
            - focus_distance * camera_z;

        Self {
            origin,
            focus_plane_horizontal,
            focus_plane_vertical,
            focus_plane_lower_left_corner,
            // SE TODO: do we have to store these, instead of normalizing horz and vert?
            camera_x,
            camera_y,
            lens_radius: aperture / 2.0,
        }
    }

    pub fn get_ray(&self, u: f64, v: f64) -> Ray {
        let jitter = self.lens_radius * Vec3::random_in_unit_disc();
        let position_on_lens =
            self.origin + self.camera_x * jitter.x() + self.camera_y * jitter.y();
        let position_on_focus_plane = self.focus_plane_lower_left_corner
            + u * self.focus_plane_horizontal
            + v * self.focus_plane_vertical;

        Ray::new(position_on_lens, position_on_focus_plane - position_on_lens)
    }
}
