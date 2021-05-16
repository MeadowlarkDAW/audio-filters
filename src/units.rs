use core::f32;
use core::f64;

use num_complex::Complex;

use libm::*;

/// Used to implement conversions to the Hertz struct
pub trait Units<T> {
    /// From hertz
    fn to_range(self, bottom: T, top: T) -> T;
    fn from_range(self, bottom: T, top: T) -> T;
    fn db_to_lin(self) -> T;
    fn lin_to_db(self) -> T;
    fn sign(self, b: T) -> T;
    fn bw_to_q(self, f0: T, fs: T) -> T;
}

impl Units<f64> for f64 {
    fn to_range(self, bottom: f64, top: f64) -> f64 {
        self * (top - bottom) + bottom
    }
    fn from_range(self, bottom: f64, top: f64) -> f64 {
        (self - bottom) / (top - bottom)
    }
    fn db_to_lin(self) -> f64 {
        pow(10.0f64, self * 0.05)
    }
    fn lin_to_db(self) -> f64 {
        log10(self.max(0.0)) * 20.0
    }
    fn sign(self, b: f64) -> f64 {
        if b < 0.0 {
            -self
        } else {
            self
        }
    }
    fn bw_to_q(self, _f0: f64, _fs: f64) -> f64 {
        1.0 / (2.0 * sinh(f64::consts::LN_2 / 2.0 * self))
    }
}

impl Units<f32> for f32 {
    //Just a copy of the f64 version with units swapped
    fn to_range(self, bottom: f32, top: f32) -> f32 {
        self * (top - bottom) + bottom
    }
    fn from_range(self, bottom: f32, top: f32) -> f32 {
        (self - bottom) / (top - bottom)
    }
    fn db_to_lin(self) -> f32 {
        powf(10.0f32, self * 0.05)
    }
    fn lin_to_db(self) -> f32 {
        log10f(self.max(0.0)) * 20.0
    }
    fn sign(self, b: f32) -> f32 {
        if b < 0.0 {
            -self
        } else {
            self
        }
    }
    fn bw_to_q(self, _f0: f32, _fs: f32) -> f32 {
        1.0 / (2.0 * sinhf(f32::consts::LN_2 / 2.0 * self))
    }
}

#[derive(Copy, Clone, Debug)]
pub struct ZSample<T> {
    pub z: Complex<T>,
    pub pow1: Complex<T>,
    pub pow2: Complex<T>,
}

impl ZSample<f64> {
    pub fn new(f_hz: f64, fs: f64) -> ZSample<f64> {
        let z = -f64::consts::TAU * f_hz / fs;
        let z = cos(z) + sin(z) * Complex::new(0.0, 1.0);
        ZSample {
            z,
            pow1: z,
            pow2: z * z,
        }
    }
}

impl ZSample<f32> {
    pub fn new(f_hz: f32, fs: f32) -> ZSample<f32> {
        let z = -f32::consts::TAU * f_hz / fs;
        let z = cosf(z) + sinf(z) * Complex::new(0.0, 1.0);
        ZSample {
            z,
            pow1: z,
            pow2: z * z,
        }
    }
}

pub fn butterworth_cascade_q(filter_order: u32, pole: u32) -> f64 {
    let mut pole = pole;
    let pole_inc = f64::consts::PI / (filter_order as f64);
    let even_order = filter_order % 2 == 0;

    let first_angle = if even_order {
        pole_inc * 0.5
    } else {
        if pole == 0 {
            return 0.5; //Also needs to be 1 pole (not biquad)
        }
        pole -= 1;
        pole_inc
    };

    1.0 / (2.0 * cos(first_angle + pole as f64 * pole_inc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_butterworth_cascade_q() {
        assert_eq!(0.7071067811865475, butterworth_cascade_q(2, 0));

        assert_eq!(0.5, butterworth_cascade_q(3, 0));
        assert_eq!(0.9999999999999998, butterworth_cascade_q(3, 1));

        assert_eq!(0.541196100146197, butterworth_cascade_q(4, 0));
        assert_eq!(1.3065629648763764, butterworth_cascade_q(4, 1));

        assert_eq!(0.5, butterworth_cascade_q(5, 0));
        assert_eq!(0.6180339887498948, butterworth_cascade_q(5, 1));
        assert_eq!(1.6180339887498947, butterworth_cascade_q(5, 2));

        assert_eq!(0.5176380902050415, butterworth_cascade_q(6, 0));
        assert_eq!(0.7071067811865475, butterworth_cascade_q(6, 1));
        assert_eq!(1.931851652578135, butterworth_cascade_q(6, 2));
    }
}
