use core::ops::{Add, Mul, Sub};
use num_complex::Complex;
use num_traits::{Float, FloatConst, One, Zero};

use crate::{
    filter_band::{BandType, FilterBand},
    units::ZSample,
};

#[derive(Copy, Clone, Debug)]
pub struct StereoFilterBand<T> {
    left: FilterBand<T>,
    right: FilterBand<T>,
}

impl<T> StereoFilterBand<T>
where
    T: Float,
    T: Zero,
    T: One,
    T: FloatConst,
    f32: Into<T>,
    u8: Into<T>,
    T: Add<Complex<T>, Output = Complex<T>>,
    T: Mul<Complex<T>, Output = Complex<T>>,
    T: Sub<Complex<T>, Output = Complex<T>>,
{
    pub fn new(sample_rate: T) -> Self {
        StereoFilterBand {
            left: FilterBand::new(sample_rate),
            right: FilterBand::new(sample_rate),
        }
    }

    pub fn update(
        &mut self,
        kind: BandType,
        in_freq: T,
        in_gain: T,
        in_bw_value: T,
        slope: T,
        sample_rate: T,
    ) {
        self.left
            .update(kind, in_freq, in_gain, in_bw_value, slope, sample_rate);
        self.right.mimic_band(&self.left);
    }

    pub fn process(&mut self, l: T, r: T) -> [T; 2] {
        [self.left.process(l), self.right.process(r)]
    }

    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        self.left.get_bode_sample(z)
    }
}
