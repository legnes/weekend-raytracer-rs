use rand::Rng;
use std::f64::consts::PI;
use std::fmt;
use std::fmt::Display;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Range, Sub, SubAssign,
};

#[derive(Clone, Copy)]
pub struct Vec3 {
    e: [f64; 3],
}

pub type Point3 = Vec3;
pub type Color = Vec3;

impl Vec3 {
    pub fn new(e0: f64, e1: f64, e2: f64) -> Self {
        Self { e: [e0, e1, e2] }
    }

    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    pub fn one() -> Self {
        Self::new(1.0, 1.0, 1.0)
    }

    pub fn x(self) -> f64 {
        self[0]
    }

    pub fn y(self) -> f64 {
        self[1]
    }

    pub fn z(self) -> f64 {
        self[2]
    }

    pub fn dot(self, other: Self) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }

    pub fn length(self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn cross(self, other: Self) -> Self {
        Self {
            e: [
                self[1] * other[2] - self[2] * other[1],
                self[2] * other[0] - self[0] * other[2],
                self[0] * other[1] - self[1] * other[0],
            ],
        }
    }

    pub fn normalized(self) -> Self {
        self / self.length()
    }

    pub fn format_color(self, divisor: u64) -> String {
        // First divide by number of samples
        let r = ((self[0] / (divisor as f64))
            // Gamma correction
            .powf(1.0 / 2.0)
            // Clamp
            .clamp(0.0, 0.999)
            // Map to byte
            * 256.0) as u64;

        let g = ((self[1] / (divisor as f64))
            .powf(1.0 / 2.0)
            .clamp(0.0, 0.999)
            * 256.0) as u64;

        let b = ((self[2] / (divisor as f64))
            .powf(1.0 / 2.0)
            .clamp(0.0, 0.999)
            * 256.0) as u64;

        format!("{} {} {}", r, g, b)
    }

    pub fn from_spherical(r: f64, phi: f64, theta: f64) -> Self {
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let cos_theta = theta.cos();
        let sin_theta = theta.sin();

        Self {
            e: [
                r * sin_theta * cos_phi,
                r * sin_theta * sin_phi,
                r * cos_theta,
            ],
        }
    }

    pub fn random() -> Self {
        Self::random_range(0.0..1.0)
    }

    pub fn random_range(range: Range<f64>) -> Self {
        let mut rng = rand::thread_rng();

        Self {
            e: [
                rng.gen_range(range.clone()),
                rng.gen_range(range.clone()),
                rng.gen_range(range.clone()),
            ],
        }
    }

    pub fn random_on_unit_sphere() -> Self {
        let mut rng = rand::thread_rng();

        let phi = rng.gen_range(0.0..(2.0 * PI));
        // Polar angle calc from https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
        let theta = f64::acos(rng.gen_range(-1.0..1.0));

        Self::from_spherical(1.0, phi, theta)
    }

    pub fn random_in_unit_sphere() -> Self {
        let mut rng = rand::thread_rng();

        let r = f64::powf(rng.gen(), 1.0 / 3.0);
        let phi = rng.gen_range(0.0..(2.0 * PI));
        let theta = f64::acos(rng.gen_range(-1.0..1.0));

        Self::from_spherical(r, phi, theta)
    }

    pub fn random_in_unit_disc() -> Self {
        let mut rng = rand::thread_rng();

        let r = rng.gen::<f64>().sqrt();
        let theta = rng.gen::<f64>() * 2.0 * PI;

        Self::new(r * theta.cos(), r * theta.sin(), 0.0)
    }

    pub fn near_zero(self) -> bool {
        const EPS: f64 = 1.0e-8;
        self[0].abs() < EPS && self[1].abs() < EPS && self[2].abs() < EPS
    }

    // Reflect across a unit normal vector n that points against self
    // Calculate the proj of self onto n in the direction of n and add it twice
    // Self is the incident vector
    pub fn reflect(self, n: Self) -> Self {
        self - 2.0 * self.dot(n) * n
    }

    // Derived from vector form of snell's law
    // Good explanation here:
    // http://cosinekitty.com/raytrace/chapter09_refraction.html
    // Self is the incident vector
    pub fn refract(self, n: Self, eta_i_over_eta_t: f64) -> Self {
        // SE TODO: min is a shortcut to normalized?
        let cos_theta = (-self).dot(n).min(1.0);
        let refracted_perpendicular = eta_i_over_eta_t * (self + cos_theta * n);
        // SE TODO: why abs? Is this picking the correct root?
        let refracted_parallel = -(1.0 - refracted_perpendicular.length().powi(2))
            .abs()
            .sqrt()
            * n;
        refracted_perpendicular + refracted_parallel
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 {
        &self.e[index]
    }
}

impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        &mut self.e[index]
    }
}

impl Add<Vec3> for Vec3 {
    // Note: changed this an other refs to Vec3 --> Self, since that's what docs do
    // Also changed "other" --> "rhs" as in docs
    type Output = Self;

    // SE TODO: Why does this have ownership of self?
    // Answer: It implements Copy, so is stored on the stack and passed by copy
    fn add(self, rhs: Self) -> Self {
        Self {
            e: [self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2]],
        }
    }
}

impl AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, rhs: Self) -> () {
        *self = Self {
            e: [self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2]],
        }
    }
}

// NOTE: Added this f64 version for me
impl Add<f64> for Vec3 {
    type Output = Self;

    fn add(self, rhs: f64) -> Self {
        Self {
            e: [self[0] + rhs, self[1] + rhs, self[2] + rhs],
        }
    }
}

impl AddAssign<f64> for Vec3 {
    fn add_assign(&mut self, rhs: f64) -> () {
        *self = Self {
            e: [self[0] + rhs, self[1] + rhs, self[2] + rhs],
        }
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            e: [self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2]],
        }
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) -> () {
        *self = Self {
            e: [self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2]],
        }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            e: [self[0] * rhs, self[1] * rhs, self[2] * rhs],
        }
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, rhs: f64) {
        *self = Self {
            e: [self[0] * rhs, self[1] * rhs, self[2] * rhs],
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        Vec3 {
            e: [self * rhs[0], self * rhs[1], self * rhs[2]],
        }
    }
}

impl Mul for Vec3 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            e: [self[0] * rhs[0], self[1] * rhs[1], self[2] * rhs[2]],
        }
    }
}

impl MulAssign for Vec3 {
    fn mul_assign(&mut self, rhs: Self) -> () {
        *self = Self {
            e: [self[0] * rhs[0], self[1] * rhs[1], self[2] * rhs[2]],
        }
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            e: [self[0] / rhs, self[1] / rhs, self[2] / rhs],
        }
    }
}

impl DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, rhs: f64) {
        *self = Self {
            e: [self[0] / rhs, self[1] / rhs, self[2] / rhs],
        }
    }
}

impl Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self {
        (-1.0) * self
    }
}

impl Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} {} {})", self[0], self[1], self[2])
    }
}
