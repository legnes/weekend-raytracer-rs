mod camera;
mod hit;
mod material;
mod ray;
mod sphere;
mod vec;

use camera::Camera;
use hit::{Hit, World};
#[allow(unused_imports)]
use material::{Dielectric, Lambertian, Metal, Scatter};
use rand::Rng;
use ray::Ray;
use sphere::Sphere;
use std::io::{stderr, Write};
use std::rc::Rc;
use std::time::Instant;
use vec::{Color, Point3};

use crate::vec::Vec3;

// SE TODO: Should try profiling this...seems to slow down...?

fn ray_color(ray: &Ray, world: &World, depth: u64) -> Color {
    if depth <= 0 {
        // Too many bounces! Assume all energy lost
        return Color::zero();
    }

    // t_min prevents hitting very near surfaces, aka shadow acne
    if let Some(hit) = world.hit(ray, 0.001, f64::INFINITY) {
        if let Some((attenuation, reflected)) = hit.material.scatter(ray, &hit) {
            attenuation * ray_color(&reflected, world, depth - 1)
        } else {
            Color::zero()
        }
    } else {
        // Background color
        let unit_direction = ray.direction().normalized();
        let t = (unit_direction.y() + 1.0) * 0.5;
        // Radial gradient around the vertical
        let from = Color::new(1.0, 1.0, 1.0);
        let to = Color::new(0.5, 0.7, 1.0);
        (1.0 - t) * from + t * to
    }
}

#[allow(dead_code)]
fn small_scene() -> World {
    let mut world = World::new();

    let mat_ground = Rc::new(Lambertian::new(Color::new(0.8, 0.8, 0.0)));
    let mat_center = Rc::new(Lambertian::new(Color::new(0.1, 0.2, 0.5)));
    let mat_left = Rc::new(Dielectric::new(1.5));
    let mat_left_inner = Rc::new(Dielectric::new(1.5));
    let mat_right = Rc::new(Metal::new(Color::new(0.8, 0.6, 0.2), 0.0));

    let sphere_ground = Sphere::new(Point3::new(0.0, -100.5, -1.0), 100.0, mat_ground);
    let sphere_center = Sphere::new(Point3::new(0.0, 0.0, -1.0), 0.5, mat_center);
    let sphere_left = Sphere::new(Point3::new(-1.0, 0.0, -1.0), 0.5, mat_left);
    let sphere_left_inner = Sphere::new(Point3::new(-1.0, 0.0, -1.0), -0.45, mat_left_inner);
    let sphere_right = Sphere::new(Point3::new(1.0, 0.0, -1.0), 0.5, mat_right);

    world.push(Box::new(sphere_ground));
    world.push(Box::new(sphere_center));
    world.push(Box::new(sphere_left));
    world.push(Box::new(sphere_left_inner));
    world.push(Box::new(sphere_right));

    world
}

#[allow(dead_code)]
fn random_scene() -> World {
    let mut rng = rand::thread_rng();
    let mut world = World::new();

    let ground_mat = Rc::new(Lambertian::new(Color::new(0.5, 0.5, 0.5)));
    let ground_sphere = Sphere::new(Point3::new(0.0, -1000.0, 0.0), 1000.0, ground_mat);
    world.push(Box::new(ground_sphere));

    for x in -11..11 {
        for z in -11..11 {
            let dx = rng.gen_range(0.0..0.9);
            let dz = rng.gen_range(0.0..0.9);
            let center = Point3::new((x as f64) + dx, 0.2, (z as f64) + dz);

            let mat: Rc<dyn Scatter>;
            let mat_sample: f64 = rng.gen();
            if mat_sample < 0.8 {
                // Diffuse
                let albedo = Color::random() * Color::random();
                mat = Rc::new(Lambertian::new(albedo));
            } else if mat_sample < 0.95 {
                // Metal
                let albedo = Color::random_range(0.4..1.0);
                let fuzz = rng.gen_range(0.0..0.5);
                mat = Rc::new(Metal::new(albedo, fuzz));
            } else {
                // Dielectric
                // SE TODO: play with ref index?
                mat = Rc::new(Dielectric::new(1.5));
            }

            let sphere = Sphere::new(center, 0.2, mat);
            world.push(Box::new(sphere));
        }
    }

    let dielectric_mat = Rc::new(Dielectric::new(1.5));
    let lambertian_mat = Rc::new(Lambertian::new(Color::new(0.4, 0.2, 0.1)));
    let metal_mat = Rc::new(Metal::new(Color::new(0.7, 0.6, 0.5), 0.0));

    let dielectric_sphere = Sphere::new(Point3::new(0.0, 1.0, 0.0), 1.0, dielectric_mat);
    let lambertian_sphere = Sphere::new(Point3::new(-4.0, 1.0, 0.0), 1.0, lambertian_mat);
    let metal_sphere = Sphere::new(Point3::new(4.0, 1.0, 0.0), 1.0, metal_mat);

    world.push(Box::new(dielectric_sphere));
    world.push(Box::new(lambertian_sphere));
    world.push(Box::new(metal_sphere));

    world
}

fn main() {
    let start = Instant::now();

    // Image
    const ASPECT_RATIO: f64 = 3.0 / 2.0;
    const IMAGE_WIDTH: u64 = 1200;
    const IMAGE_HEIGHT: u64 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u64;
    const SAMPLES_PER_PIXEL: u64 = 500;
    const MAX_DEPTH: u64 = 50;

    // World
    let world = random_scene();

    // Camera
    let look_from = Point3::new(13.0, 2.0, 3.0);
    let look_at = Point3::new(0.0, 0.0, 0.0);
    let world_up = Vec3::new(0.0, 1.0, 0.0);
    let focus_distance = 10.0;
    let fov_degrees = 20.0;
    let aperture = 0.1;

    let camera = Camera::new(
        look_from,
        look_at,
        world_up,
        fov_degrees,
        ASPECT_RATIO,
        aperture,
        focus_distance,
    );

    // Output Format
    // ASCII
    println!("P3");
    // Dimensions
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    // Max color
    println!("255");

    // Pixel values
    let mut rng = rand::thread_rng();
    for j in (0..IMAGE_HEIGHT).rev() {
        eprint!(
            "\rRunning scanline: {:3} of {}",
            IMAGE_HEIGHT - j,
            IMAGE_HEIGHT
        );
        stderr().flush().unwrap();

        for i in 0..IMAGE_WIDTH {
            let mut pixel_color = Color::zero();
            for _ in 0..SAMPLES_PER_PIXEL {
                // SE TODO: Try adding stratification
                let random_u: f64 = rng.gen();
                let random_v: f64 = rng.gen();

                let u = (i as f64 + random_u) / ((IMAGE_WIDTH - 1) as f64);
                let v = (j as f64 + random_v) / ((IMAGE_HEIGHT - 1) as f64);

                let ray = camera.get_ray(u, v);
                pixel_color += ray_color(&ray, &world, MAX_DEPTH)
            }

            println!("{}", pixel_color.format_color(SAMPLES_PER_PIXEL));
        }
    }
    eprintln!("\nDone in {} seconds!", start.elapsed().as_secs());
}
