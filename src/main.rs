mod camera;
mod hit;
mod ray;
mod sphere;
mod vec;

use camera::Camera;
use hit::{Hit, World};
use rand::Rng;
use ray::Ray;
use sphere::Sphere;
use std::io::{stderr, Write};
use vec::{Color, Point3};

// SE TODO: You made it to "Antialiasing"

fn ray_color(ray: &Ray, world: &World) -> Color {
    if let Some(hit) = world.hit(ray, 0.0, f64::INFINITY) {
        hit.normal * 0.5 + 0.5
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
    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: u64 = 256;
    const IMAGE_HEIGHT: u64 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u64;
    const SAMPLES_PER_PIXEL: u64 = 100;

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
                pixel_color += ray_color(&ray, &world)
            }

            println!("{}", pixel_color.format_color(SAMPLES_PER_PIXEL));
        }
    }
    eprintln!("\nDone!");
}
