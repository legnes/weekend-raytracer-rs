use rand::Rng;
use std::f64::consts::PI;
use std::fmt;
use std::fmt::Display;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Range, Sub, SubAssign,
};

#[derive(Clone, Copy, Debug)]
pub struct Vec3 {
    e: [f64; 3],
}

pub type Point3 = Vec3;
pub type Color = Vec3;

impl Color {
    pub fn debug() -> Self {
        Self::new(1.0, 0.0, 1.0)
    }

    pub fn display(self) -> (u8, u8, u8) {
        // Gamma correction, then clamp, then map to byte
        let r = (self[0].powf(1.0 / 2.2).clamp(0.0, 0.999) * 256.0) as u8;
        let g = (self[1].powf(1.0 / 2.2).clamp(0.0, 0.999) * 256.0) as u8;
        let b = (self[2].powf(1.0 / 2.2).clamp(0.0, 0.999) * 256.0) as u8;

        (r, g, b)
    }

    pub fn format(self) -> String {
        let (r, g, b) = self.display();

        format!("{} {} {}", r, g, b)
    }

    pub fn luminance(self) -> f64 {
        // https://en.wikipedia.org/wiki/Relative_luminance
        self.dot(Vec3::new(0.2126, 0.7152, 0.0722))

        // https://www.w3.org/TR/AERT/#color-contrast
        // self.dot(Vec3::new(0.299, 0.587, 0.114))
    }

    pub fn expose(self, d_stop: f64) -> Self {
        self * f64::powf(2.0, d_stop)
    }

    // http://filmicworlds.com/blog/filmic-tonemapping-operators/
    fn tonemap_uncharted_2_helper(self) -> Self {
        let a = 0.15;
        let b = 0.50;
        let c = 0.10;
        let d = 0.20;
        let e = 0.02;
        let f = 0.30;
        ((self * (self * a + c * b) + d * e) / (self * (self * a + b) + d * f)) - e / f
    }

    pub fn tonemap_uncharted_2(self) -> Self {
        let exposure_bias = 2.0;
        let exposed = self * exposure_bias;
        let curr = exposed.tonemap_uncharted_2_helper();

        let white_scale = Vec3::one() / Vec3::fill(11.2).tonemap_uncharted_2_helper();
        curr * white_scale
    }

    pub fn tonemap_hejl_burgess_dawson(self) -> Self {
        let x = (self - 0.004).max(Vec3::zero());
        (x * (x * 6.2 + 0.5)) / (x *(x * 6.2 + 1.7) + 0.06)
    }

    // https://64.github.io/tonemapping/
    pub fn tonemap_clamp(self) -> Self {
        self.clamp(0.0, 1.0)
    }

    pub fn tonemap_reinhard(self) -> Self {
        self / (self + 1.0)
    }

    pub fn tonemap_reinhard_luminance(self) -> Self {
        self / (self.luminance() + 1.0)
    }

    pub fn tonemap_aces_approximate(self) -> Self {
        let temp = self * 0.6;
        let a = 2.51;
        let b = 0.03;
        let c = 2.43;
        let d = 0.59;
        let e = 0.14;
        (temp * (temp * a + b)) / (temp * (temp * c + d) + e).clamp(0.0, 1.0)
    }
}

impl Vec3 {
    pub fn new(e0: f64, e1: f64, e2: f64) -> Self {
        Self { e: [e0, e1, e2] }
    }

    pub fn fill(v: f64) -> Self {
        Self::new(v, v, v)
    }

    pub fn zero() -> Self {
        Self::fill(0.0)
    }

    pub fn one() -> Self {
        Self::fill(1.0)
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

    pub fn min(self, other: Self) -> Self {
        Self {
            e: [
                self[0].min(other[0]),
                self[1].min(other[1]),
                self[2].min(other[2]),
            ],
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            e: [
                self[0].max(other[0]),
                self[1].max(other[1]),
                self[2].max(other[2]),
            ],
        }
    }

    pub fn clamp(self, min: f64, max: f64) -> Self {
        Self {
            e: [
                self[0].clamp(min, max),
                self[1].clamp(min, max),
                self[2].clamp(min, max),
            ],
        }
    }

    pub fn is_nan(self) -> bool {
        f64::is_nan(self.e[0]) || f64::is_nan(self.e[1]) || f64::is_nan(self.e[2])
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

// NOTE: Added this f64 version for me
impl Sub<f64> for Vec3 {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self {
        Self {
            e: [self[0] - rhs, self[1] - rhs, self[2] - rhs],
        }
    }
}

impl SubAssign<f64> for Vec3 {
    fn sub_assign(&mut self, rhs: f64) -> () {
        *self = Self {
            e: [self[0] - rhs, self[1] - rhs, self[2] - rhs],
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

impl Div for Vec3 {
    type Output = Self;

    fn div(self, rhs: Vec3) -> Self::Output {
        Self {
            e: [self[0] / rhs[0], self[1] / rhs[1], self[2] / rhs[2]],
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

impl PartialEq for Vec3 {
    fn eq(&self, other: &Self) -> bool {
        self.e[0] == other.e[0] && self.e[1] == other.e[1] && self.e[2] == other.e[2]
    }
}

impl PartialEq<f64> for Vec3 {
    fn eq(&self, other: &f64) -> bool {
        self.e[0] == *other && self.e[1] == *other && self.e[2] == *other
    }
}
