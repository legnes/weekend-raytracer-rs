#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;
// https://docs.rs/dhat/latest/dhat/
// https://nnethercote.github.io/dh_view/dh_view.html

// Credits
// https://misterdanb.github.io/raytracinginrust/
// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/

mod bvh;
mod camera;
mod hit;
mod material;
mod ray;
mod sphere;
mod vec;

use camera::Camera;
use hit::{Hit, Scene, World};
#[allow(unused_imports)]
use material::{Dielectric, Lambertian, Metal, Scatter};
use rand::Rng;
use ray::Ray;
#[allow(unused_imports)]
use rayon::prelude::*;
use sphere::Sphere;
use std::io::{stderr, Write};
use std::sync::Arc;
use std::time::Instant;
use vec::{Color, Point3, Vec3};

// SE TODO: Next steps
//
//  - Make bvh faster (keep following tutorial)
//  - Rainbows - Make dielectric scatter frequency-dependence and shoot 1 frequency per ray
//  - HDR rendering
//  - PBR materials
//  - Image-based lighting
//  - WebGPU version
//  - Wasm version
//
//  - Lights — You can do this explicitly, by sending shadow rays to lights, or it can be done implicitly by making some objects emit light, biasing scattered rays toward them, and then downweighting those rays to cancel out the bias. Both work. I am in the minority in favoring the latter approach.
//  - Triangles — Most cool models are in triangle form. The model I/O is the worst and almost everybody tries to get somebody else’s code to do this.
//  - Surface Textures — This lets you paste images on like wall paper. Pretty easy and a good thing to do.
//  - Solid textures — Ken Perlin has his code online. Andrew Kensler has some very cool info at his blog.
//  - Volumes and Media — Cool stuff and will challenge your software architecture. I favor making volumes have the hittable interface and probabilistically have intersections based on density. Your rendering code doesn’t even have to know it has volumes with that method.
//  - Parallelism — Run N copies of your code on N cores with different random seeds. Average the N runs. This averaging can also be done hierarchically where N/2 pairs can be averaged to get N/4 images, and pairs of those can be averaged. That method of parallelism should extend well into the thousands of cores with very little coding.
//

// SE TDOO: Done
//  - Profiling - seems to slow down...is there a leak?
//      - no obvious leak, but spend all time in sphere intersection
//  - Accelerate hittest - rtree or something?
//      - added bvh, which sped things up by an factor of 7 for 500 objects and a factor of 18 for 1800 objects
//

fn ray_color(ray: &Ray, scene: &Scene, depth: u64) -> Color {
    if depth <= 0 {
        // Too many bounces! Assume all energy lost
        return Color::zero();
    }

    // t_min prevents hitting very near surfaces, aka shadow acne
    if let Some(hit) = scene.hit(ray, 0.001, f64::INFINITY) {
        if let Some((attenuation, reflected)) = hit.material.scatter(ray, &hit) {
            attenuation * ray_color(&reflected, scene, depth - 1)
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
fn small_world() -> World {
    let mut world = World::new();

    let mat_ground = Arc::new(Lambertian::new(Color::new(0.8, 0.8, 0.0)));
    let mat_center = Arc::new(Lambertian::new(Color::new(0.1, 0.2, 0.5)));
    let mat_left = Arc::new(Dielectric::new(1.5));
    let mat_left_inner = Arc::new(Dielectric::new(1.5));
    let mat_right = Arc::new(Metal::new(Color::new(0.8, 0.6, 0.2), 0.0));

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
fn random_world(n: i32) -> World {
    let mut rng = rand::thread_rng();
    let mut world = World::new();

    let ground_mat = Arc::new(Lambertian::new(Color::new(0.5, 0.5, 0.5)));
    let ground_sphere = Sphere::new(Point3::new(0.0, -1000.0, 0.0), 1000.0, ground_mat);
    world.push(Box::new(ground_sphere));

    for x in -n..n {
        for z in -n..n {
            let dx = rng.gen_range(0.0..0.9);
            let dz = rng.gen_range(0.0..0.9);
            let center = Point3::new((x as f64) + dx, 0.2, (z as f64) + dz);

            let mat: Arc<dyn Scatter>;
            let mat_sample: f64 = rng.gen();
            if mat_sample < 0.8 {
                // Diffuse
                let albedo = Color::random() * Color::random();
                mat = Arc::new(Lambertian::new(albedo));
            } else if mat_sample < 0.95 {
                // Metal
                let albedo = Color::random_range(0.4..1.0);
                let fuzz = rng.gen_range(0.0..0.5);
                mat = Arc::new(Metal::new(albedo, fuzz));
            } else {
                // Dielectric
                // SE TODO: play with ref index?
                mat = Arc::new(Dielectric::new(1.5));
            }

            let sphere = Sphere::new(center, 0.2, mat);
            world.push(Box::new(sphere));
        }
    }

    let dielectric_mat = Arc::new(Dielectric::new(1.5));
    let lambertian_mat = Arc::new(Lambertian::new(Color::new(0.4, 0.2, 0.1)));
    let metal_mat = Arc::new(Metal::new(Color::new(0.7, 0.6, 0.5), 0.0));

    let dielectric_sphere = Sphere::new(Point3::new(0.0, 1.0, 0.0), 1.0, dielectric_mat);
    let lambertian_sphere = Sphere::new(Point3::new(-4.0, 1.0, 0.0), 1.0, lambertian_mat);
    let metal_sphere = Sphere::new(Point3::new(4.0, 1.0, 0.0), 1.0, metal_mat);

    world.push(Box::new(dielectric_sphere));
    world.push(Box::new(lambertian_sphere));
    world.push(Box::new(metal_sphere));

    world
}

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let start = Instant::now();

    // Image
    const ASPECT_RATIO: f64 = 3.0 / 2.0;
    const IMAGE_WIDTH: u64 = 256; // 1024;
    const IMAGE_HEIGHT: u64 = ((IMAGE_WIDTH as f64) / ASPECT_RATIO) as u64;
    const SAMPLES_PER_PIXEL: u64 = 100; // 500;
    const MAX_DEPTH: u64 = 10; // 50;

    // World
    let scene = Scene::new(random_world(11));

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
    // SE TODO: Run these in 4x4 tiles or something, for better data locality?
    for j in (0..IMAGE_HEIGHT).rev() {
        eprint!(
            "\rRunning scanline: {:3} of {}",
            IMAGE_HEIGHT - j,
            IMAGE_HEIGHT
        );
        stderr().flush().unwrap();

        // let scanline: Vec<Color> = (0..IMAGE_WIDTH)
        //     .into_iter()
        //     .map(|i| {
        let scanline: Vec<Color> = (0..IMAGE_WIDTH)
            .into_par_iter()
            .map(|i| {
                let mut rng = rand::thread_rng();
                let mut pixel_color = Color::zero();
                for _ in 0..SAMPLES_PER_PIXEL {
                    // SE TODO: Try adding stratification
                    let subpixel_u: f64 = rng.gen();
                    let subpixel_v: f64 = rng.gen();

                    let u = (i as f64 + subpixel_u) / ((IMAGE_WIDTH - 1) as f64);
                    let v = (j as f64 + subpixel_v) / ((IMAGE_HEIGHT - 1) as f64);

                    let ray = camera.get_ray(u, v);
                    pixel_color += ray_color(&ray, &scene, MAX_DEPTH)
                }

                pixel_color
            })
            .collect();

        for pixel_color in scanline {
            println!("{}", pixel_color.format_color(SAMPLES_PER_PIXEL));
        }
    }
    eprintln!("\nDone in {} seconds!", start.elapsed().as_secs());
}
