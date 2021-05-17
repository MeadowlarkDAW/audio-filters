use core::ops::{Add, Mul, Sub};

use num_complex::Complex;

use num_traits::{Float, FloatConst, NumCast, One, Zero};

pub trait FP:
    Sized
    + Copy
    + Float
    + Zero
    + One
    + FloatConst
    + From<f32>
    + From<u8>
    + Into<Complex<Self>>
    + Add<Complex<Self>, Output = Complex<Self>>
    + Mul<Complex<Self>, Output = Complex<Self>>
    + Sub<Complex<Self>, Output = Complex<Self>>
{
}

impl<T> FP for T
where
    T: Sized,
    T: Copy,
    T: Float,
    T: Zero,
    T: One,
    T: FloatConst,
    T: From<f32>,
    T: From<u8>,
    T: Into<Complex<Self>>,
    T: Add<Complex<Self>, Output = Complex<Self>>,
    T: Mul<Complex<Self>, Output = Complex<Self>>,
    T: Sub<Complex<Self>, Output = Complex<Self>>,
{
}

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

impl<T: FP> Units<T> for T {
    fn to_range(self, bottom: T, top: T) -> T {
        self * (top - bottom) + bottom
    }
    fn from_range(self, bottom: T, top: T) -> T {
        (self - bottom) / (top - bottom)
    }
    fn db_to_lin(self) -> T {
        Into::<T>::into(10.0).powf(self * Into::<T>::into(0.05))
    }
    fn lin_to_db(self) -> T {
        (self.max(T::zero())).log10() * Into::<T>::into(20.0)
    }
    fn sign(self, b: T) -> T {
        if b < T::zero() {
            -self
        } else {
            self
        }
    }
    fn bw_to_q(self, _f0: T, _fs: T) -> T {
        let two = NumCast::from(2.0).unwrap();
        T::one() / (two * (T::LN_2() / two * self).sinh())
    }
}

#[derive(Copy, Clone, Debug)]
pub struct ZSample<T> {
    pub z: Complex<T>,
    pub pow1: Complex<T>,
    pub pow2: Complex<T>,
}

impl<T: FP> ZSample<T> {
    pub fn new(f_hz: T, fs: T) -> ZSample<T> {
        let z = -T::TAU() * f_hz / fs;
        let z: Complex<T> =
            z.cos().into() + z.sin().into() * Complex::<T>::new(T::zero(), T::one());
        ZSample {
            z,
            pow1: z,
            pow2: z * z,
        }
    }
}

pub fn butterworth_cascade_q<T: FP>(filter_order: u8, pole: u8) -> T {
    let mut pole = pole;
    let pole_inc: T = T::PI() / (NumCast::from(filter_order).unwrap());
    let even_order = filter_order & 1 == 0;
    let point_five = NumCast::from(0.5).unwrap();
    let two: T = NumCast::from(2.0).unwrap();

    let first_angle = if even_order {
        pole_inc * point_five
    } else {
        if pole == 0 {
            return point_five; //Also needs to be 1 pole (not biquad)
        }
        pole -= 1;
        pole_inc
    };
    let fpole: T = NumCast::from(pole).unwrap();
    let a: T = first_angle + fpole * pole_inc;
    T::one() / (two * a.cos())
}

pub fn test_func<T: FP>(x: T) -> T {
    x
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
        dbg!(butterworth_cascade_q::<f64>(5, 2));
    }
}
