mod camera;
mod hit;
mod ray;
mod sphere;
mod vec;
mod shading;

use camera::Camera;
use hit::{Hit, World};
use rand::Rng;
use ray::Ray;
use shading::DiffuseModel;
use sphere::Sphere;
use std::io::{stderr, Write};
use std::time::Instant;
use vec::{Color, Point3};

// SE TODO: You made it to "True Lambertian Reflection"
// SE TODO: Something seems wrong...

fn ray_color(ray: &Ray, world: &World, depth: u64) -> Color {
    if depth <= 0 {
        // Too many bounces! Assume all energy lost
        return Color::zero();
    }

    // t_min prevents hitting very near surfaces, aka shadow acne
    if let Some(hit) = world.hit(ray, 0.001, f64::INFINITY) {
        let reflection = shading::diffuse(hit, DiffuseModel::Lambert);
        // Decrease energy for each bounce.
        0.5 * ray_color(&reflection, world, depth - 1)
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

fn main() {
    let start = Instant::now();

    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u64 = 256;
    const IMAGE_HEIGHT: u64 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u64;
    const SAMPLES_PER_PIXEL: u64 = 100;
    const MAX_DEPTH: u64 = 5;

    // World
    let mut world = World::new();
    world.push(Box::new(Sphere::new(Point3::new(0.0, 0.0, -1.0), 0.5)));
    world.push(Box::new(Sphere::new(Point3::new(0.0, -100.5, -1.0), 100.0)));

    // Camera
    let camera = Camera::new();

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
        // SE TODO: had to change this line...?
        eprint!("\rRunning scanline: {:3}", IMAGE_HEIGHT - j + 1);
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
