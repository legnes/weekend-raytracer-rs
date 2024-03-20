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
mod light;
mod material;
mod ray;
mod sphere;
mod vec;

use camera::Camera;
use hit::{Hit, Scene, World};
use light::{DirectionalLight, Light};
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

// For now we are operating in spectral radiance
const RADIANCE_SKY: f64 = 500.0 / 3.0;
const RADIANCE_SUN: f64 = 1000.0 / 3.0;

// SE TODO: Next steps
//  - Implement a few BRDFs?
//  - Track more spectral bands, make brdf wavelength dependent?
//  - Move back to hittable light sources?
//  - Add analytic sky model
//
//  - Volumes!
//  - Subsurface scattering
//  - Rainbows - Make dielectric scatter frequency-dependence and shoot 1 frequency per ray
//      - Each ray gets a wavelength
//      - Scatter/brdf depends on wavelength
//      - For shading, convert -> XYZ (mul by CIE color matching) -> rgb
//  - HDR rendering
//  - PBR materials
//  - Image-based lighting
//  - WebGPU version
//  - Wasm version
//
//  - Lights — You can do this explicitly, by sending shadow rays to lights, or it can be done implicitly by making some objects emit light, biasing scattered rays toward them, and then downweighting those rays to cancel out the bias. Both work. I am in the minority in favoring the latter approach.
//      - For more on shadow rays: https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/ligth-and-shadows.html
//      - Let's try the other version, where we just make emissive hittables in the world. How do you weight samples towards the lights? And how do you calc/cancel that weight? One option is to pick a light with some prob and sample directly towards it? Essentially shadow ray...
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
//  - Make bvh faster (keep following tutorial)
//      - added surface area heuristic to bvh construction, seemed to slow traversal down a lot (possibly bug in sah impl?)
//      - added sorted traversal (hit near child first), also seemed to slow things down
//      - removed both...
//  - Lighting
//      - prototyped shadow ray approach for directional lights. looked good, but introduced entirely parallel shading model and struggled with dielectric shadows
//      - converted to hdr, experimented with different tonemappers
//      - Added IBL from EXR -- struggled reconciling with sun radiance units?
//

fn ray_color(ray: &Ray, scene: &Scene, depth: u64, skybox: Arc<Vec<Vec<[f32; 4]>>>) -> Color {
    if depth <= 0 {
        // Too many bounces! Assume all energy lost
        return Color::zero();
    }

    // t_min prevents hitting very near surfaces, aka shadow acne
    if let Some(hit) = scene.hit(ray, 0.001, f64::INFINITY) {
        let ambient_light = if let Some((attenuation, reflected)) = hit.material.scatter(ray, &hit)
        {
            attenuation * ray_color(&reflected, scene, depth - 1, skybox.clone())
        } else {
            Color::zero()
        };

        let mut direct_light = Vec3::zero();
        for light in scene.lights() {
            let shadow_ray = light.get_shadow_ray(hit.position);
            let is_lit = if let Some(shadow_hit) = scene.hit(&shadow_ray, 0.001, f64::INFINITY) {
                !shadow_hit.material.casts_shadow()
            } else {
                true
            };
            if is_lit {
                direct_light += hit.material.shade(ray, &shadow_ray, &hit) * light.color();
            }
        }

        ambient_light + direct_light
    } else {
        // Background color
        let unit_direction = ray.direction().normalized();
        let theta = unit_direction.y().acos();
        let xz = (ray.direction() * Vec3::new(1.0, 0.0, 1.0)).normalized();
        let phi = xz.z().signum() * xz.x().acos();
        let u = (phi / (std::f64::consts::PI * 2.0)).clamp(-1.0, 1.0) - 1.0 * -0.5;
        let v = (theta / (std::f64::consts::PI * 1.0)).clamp(0.0, 1.0);
        let row = (v * (skybox.len() - 1) as f64).floor() as usize;
        let col = (u * (skybox[0].len() - 1) as f64).floor() as usize;
        let pixel = skybox[row][col];
        RADIANCE_SKY * Color::new(pixel[0] as f64, pixel[1] as f64, pixel[2] as f64)
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

#[allow(dead_code)]
fn random_lights() -> Vec<Box<dyn Light>> {
    let sun1 = Box::new(DirectionalLight::new(
        Vec3::new(1.0, -1.0, 1.0),
        Color::one() * RADIANCE_SUN,
    ));

    let sun2 = Box::new(DirectionalLight::new(
        Vec3::random_range(-1.0..-0.1),
        Color::random_range(0.25..1.0) * RADIANCE_SUN,
    ));

    vec![sun1, sun2]
}

#[allow(dead_code)]
fn static_lights() -> Vec<Box<dyn Light>> {
    let sun1 = Box::new(DirectionalLight::new(
        Vec3::new(1.0, -1.0, 1.0),
        Color::new(0.8, 1.0, 0.9) * RADIANCE_SUN,
    ));

    let sun2 = Box::new(DirectionalLight::new(
        Vec3::new(-0.4, -1.0, 1.0),
        Color::new(1.0, 0.9, 0.8) * RADIANCE_SUN / 2.0,
    ));

    vec![sun1, sun2]
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

    // Scene
    let scene = Scene::new(random_world(11), random_lights());

    // Skybox
    let skybox = exr::prelude::read_first_rgba_layer_from_file(
        "./assets/venice_sunset_4k.exr",
        |resolution, _| {
            let default_pixel = [0.0, 0.0, 0.0, 0.0];
            let empty_line = vec![default_pixel; resolution.width()];
            let empty_image = vec![empty_line; resolution.height()];
            empty_image
        },
        // transfer the colors from the file to your image type,
        // NOTE: it seems like openExr stores in linear space? But some parsers apply gamma correction...
        |pixel_vector, position, (r, g, b, a): (f32, f32, f32, f32)| {
            pixel_vector[position.y()][position.x()] = [r, g, b, a]
        },
    )
    .expect("should get a skybox");
    let skybox = Arc::new(skybox.layer_data.channel_data.pixels);

    // Camera
    let look_from = Point3::new(13.0, 2.0, 3.0);
    let look_at = Point3::new(0.0, 1.0, 0.0);
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
                    pixel_color += ray_color(&ray, &scene, MAX_DEPTH, skybox.clone())
                }

                pixel_color / (SAMPLES_PER_PIXEL as f64)
            })
            .collect();

        for pixel_color in scanline {
            println!(
                "{}",
                pixel_color.expose(-8.0).tonemap_aces_approximate().format()
            );
        }
    }
    eprintln!("\nDone in {} seconds!", start.elapsed().as_secs());
}
