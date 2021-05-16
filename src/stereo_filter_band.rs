use num_complex::Complex;

use crate::{
    filter_band::{BandType, FilterBand},
    units::ZSample,
};

#[derive(Copy, Clone, Debug)]
pub struct StereoFilterBand {
    left: FilterBand,
    right: FilterBand,
}

impl StereoFilterBand {
    pub fn new(sample_rate: f64) -> Self {
        StereoFilterBand {
            left: FilterBand::new(sample_rate),
            right: FilterBand::new(sample_rate),
        }
    }

    pub fn update(
        &mut self,
        kind: BandType,
        in_freq: f64,
        in_gain: f64,
        in_bw_value: f64,
        slope: f64,
        sample_rate: f64,
    ) {
        self.left
            .update(kind, in_freq, in_gain, in_bw_value, slope, sample_rate);
        self.right.mimic_band(&self.left);
    }

    pub fn process(&mut self, l: f64, r: f64) -> [f64; 2] {
        [self.left.process(l), self.right.process(r)]
    }

    pub fn get_bode_sample(&self, z: ZSample<f64>) -> Complex<f64> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        self.left.get_bode_sample(z)
    }
}
