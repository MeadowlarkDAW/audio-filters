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
#[cfg(test)]
mod tests {
    use super::*;

    fn rand(x: f32) -> f32 {
        ((x * 12.9898).sin() * 43758.5453).fract()
    }

    #[test]
    fn it_works() {
        let mut left: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();
        let mut right: Vec<f32> = (1000..2000).map(|x| rand(x as f32)).collect();

        let sample_rate = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let bandwidth = 1.0;
        let slope = 4.0;
        let mut filter = StereoFilterBand::new(sample_rate);
        filter.update(BandType::HighShelf, f0, gain, bandwidth, slope, sample_rate);
        for i in 0..1000 {
            let [l_out, r_out] = filter.process(left[i], right[i]);
            left[i] = l_out;
            right[i] = r_out;
        }
    }
}
