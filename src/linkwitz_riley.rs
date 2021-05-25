use num_complex::Complex;

use crate::{
    filter_band::{FilterBand, FilterBandCoefficients, ProcessType},
    units::{ZSample, FP},
};

#[derive(Copy, Clone, Debug)]
pub struct LinkwitzRileyCoefficients<T: FP> {
    pub coeffs: FilterBandCoefficients<T>,
}

impl<T: FP> LinkwitzRileyCoefficients<T> {
    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        self.coeffs.get_bode_sample(z) * self.coeffs.get_bode_sample(z)
    }

    //The resulting Linkwitz-Riley filter will have 2x to order of the input coefficients and 2x gain
    pub fn from(coeffs: FilterBandCoefficients<T>) -> Self {
        LinkwitzRileyCoefficients { coeffs }
    }
}

#[derive(Copy, Clone)]
pub struct LinkwitzRileyBand<T: FP> {
    pub filter1: FilterBand<T>,
    pub filter2: FilterBand<T>,
    pub process: fn(&mut Self, T) -> T,
}

impl<T: FP> LinkwitzRileyBand<T> {
    pub fn from(lw_coeffs: &LinkwitzRileyCoefficients<T>) -> LinkwitzRileyBand<T> {
        LinkwitzRileyBand {
            filter1: FilterBand::from(&lw_coeffs.coeffs),
            filter2: FilterBand::from(&lw_coeffs.coeffs),
            process: LinkwitzRileyBand::get_process(lw_coeffs.coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, input_sample: T) -> T {
        self.filter2
            .process_iir1_only(self.filter1.process_iir1_only(input_sample))
    }

    pub fn process_iir2_only(&mut self, input_sample: T) -> T {
        self.filter2
            .process_iir2_only(self.filter1.process_iir2_only(input_sample))
    }

    pub fn process_even_order_cascade(&mut self, input_sample: T) -> T {
        self.filter2
            .process_even_order_cascade(self.filter1.process_even_order_cascade(input_sample))
    }

    pub fn process_odd_order_cascade(&mut self, input_sample: T) -> T {
        self.filter2
            .process_odd_order_cascade(self.filter1.process_odd_order_cascade(input_sample))
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, T) -> T {
        match process_type {
            ProcessType::ProcessIIR1Only => LinkwitzRileyBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => LinkwitzRileyBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => LinkwitzRileyBand::process_even_order_cascade,
            ProcessType::ProcessOddOrderCascade => LinkwitzRileyBand::process_odd_order_cascade,
        }
    }

    pub fn update(&mut self, lw_coeffs: &LinkwitzRileyCoefficients<T>) {
        self.filter1.update(&lw_coeffs.coeffs);
        self.filter2.update(&lw_coeffs.coeffs);
        self.process = LinkwitzRileyBand::get_process(lw_coeffs.coeffs.process);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand(x: f32) -> f32 {
        ((x * 12.9898).sin() * 43758.5453).fract()
    }

    #[test]
    fn test_filter_band() {
        let mut left: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();
        let mut right: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let slope = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);
        let lw_coeffs = LinkwitzRileyCoefficients::from(coeffs);

        let mut filter_left = LinkwitzRileyBand::from(&lw_coeffs);
        let mut filter_right = LinkwitzRileyBand::from(&lw_coeffs);

        for i in 0..1000 {
            left[i] = (filter_left.process)(&mut filter_left, left[i]);
            right[i] = (filter_right.process)(&mut filter_right, right[i]);
        }

        dbg!(left[500], right[500]);
    }
}
