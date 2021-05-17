use crate::units::FP;
use num_complex::Complex;

use num_traits::NumCast;

use crate::{
    first_order_iir::{IIR1Coefficients, IIR1},
    second_order_iir::{IIR2Coefficients, IIR2},
    units::{butterworth_cascade_q, Units, ZSample},
    MAX_POLE_COUNT,
};

#[derive(Copy, Clone, Debug)]
pub struct IIRCoefficientsSet<T: FP> {
    pub iir1: IIR1Coefficients<T>,
    pub iir2: [IIR2Coefficients<T>; MAX_POLE_COUNT],
    pub iir1_enabled: bool,
}

impl<T: FP> IIRCoefficientsSet<T> {
    pub fn new(sample_rate: T) -> IIRCoefficientsSet<T> {
        let iir2_coeffs = IIR2Coefficients::bell(sample_rate, T::one(), 1000.0.into(), T::zero());
        let iir1_coeffs = IIR1Coefficients::<T>::lowpass(1000.0.into(), sample_rate);
        IIRCoefficientsSet {
            iir1_enabled: false,
            iir1: iir1_coeffs,
            iir2: [iir2_coeffs; MAX_POLE_COUNT],
        }
    }

    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        let mut y = if self.iir1_enabled {
            self.iir1.get_bode_sample(z.z) * self.iir2[0].get_bode_sample(z)
        } else {
            self.iir2[0].get_bode_sample(z)
        };
        for band_a in self.iir2.iter().skip(1) {
            //TODO only do needed iir2 bands skip others
            y = y * band_a.get_bode_sample(z);
        }
        y
    }
}

#[derive(Copy, Clone)]
pub struct FilterBand<T: FP> {
    iir1: IIR1<T>,
    iir2: [IIR2<T>; MAX_POLE_COUNT],
    pub coeffs: IIRCoefficientsSet<T>,
    f0: T,
    gain: T,
    bw_value: T,
    slope: T,
    u_slope: u8,
    odd_order: bool,
    start_pole: usize,
    iir2_cascade_count: usize,
    sample_rate: T,
    pub process: fn(&mut Self, T) -> T,
}

impl<T: FP> FilterBand<T> {
    pub fn new(sample_rate: T) -> FilterBand<T> {
        let coeffs = IIRCoefficientsSet::<T>::new(sample_rate);
        FilterBand {
            iir1: IIR1::<T>::new(coeffs.iir1),
            iir2: [IIR2::<T>::new(coeffs.iir2[0]); MAX_POLE_COUNT],
            coeffs,
            f0: T::zero(),
            gain: T::zero(),
            bw_value: T::zero(),
            slope: 2.0.into(),
            odd_order: true,
            u_slope: 2,
            start_pole: 0,
            iir2_cascade_count: 1,
            sample_rate: 48000.0.into(),
            process: FilterBand::process_iir1_only,
        }
    }

    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        self.coeffs.get_bode_sample(z)
    }

    pub fn process_iir1_only(&mut self, x: T) -> T {
        self.iir1.process(x)
    }

    pub fn process_iir2_only(&mut self, x: T) -> T {
        self.iir2[0].process(x)
    }

    pub fn process_even_order_cascade(&mut self, x: T) -> T {
        let mut x = x;
        for i in self.start_pole..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn process_odd_order_cascade(&mut self, x: T) -> T {
        let mut x = self.iir1.process(x);
        for i in self.start_pole..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn update_coeffs(&mut self, coeffs: IIRCoefficientsSet<T>) {
        for (filter, coeff) in self.iir2.iter_mut().zip(coeffs.iir2.iter()) {
            filter.update_coefficients(*coeff)
        }
        self.iir1.update_coefficients(coeffs.iir1);
    }

    pub fn mimic_band(&mut self, band: &FilterBand<T>) {
        self.coeffs = band.coeffs;
        self.f0 = band.f0;
        self.gain = band.gain;
        self.bw_value = band.bw_value;
        self.slope = band.slope;
        self.odd_order = band.odd_order;
        self.u_slope = band.u_slope;
        self.start_pole = band.start_pole;
        self.iir2_cascade_count = band.iir2_cascade_count;
        self.sample_rate = band.sample_rate;
        self.process = band.process;
        self.update_coeffs(self.coeffs);
    }

    fn set_params(&mut self, f0: T, gain: T, bw: T, slope: T, iir1_possible: bool, fs: T) {
        self.f0 = f0;
        self.gain = gain;
        self.bw_value = bw;
        self.slope = slope;
        self.sample_rate = fs;

        self.u_slope = NumCast::from(slope).unwrap();

        self.odd_order = self.u_slope & 1 == 1;

        self.start_pole = if iir1_possible && self.odd_order {
            1
        } else {
            0
        } as usize;
        self.iir2_cascade_count = ((self.u_slope as usize + self.start_pole) / 2) as usize;
    }

    pub fn lowpass(&mut self, f0: T, bw: T, slope: T, fs: T) {
        self.set_params(f0, T::zero(), bw, slope, true, fs);
        if self.odd_order {
            self.coeffs.iir1 = IIR1Coefficients::lowpass(f0, fs);
            if self.u_slope <= 1 {
                self.update_coeffs(self.coeffs);
                self.process = FilterBand::process_iir1_only;
                return;
            }
            self.process = FilterBand::process_odd_order_cascade;
        } else {
            self.process = FilterBand::process_even_order_cascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        for i in self.start_pole..self.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
            self.coeffs.iir2[i] = IIR2Coefficients::lowpass(f0, q_value * q_offset, fs);
        }
        self.update_coeffs(self.coeffs);
    }

    pub fn highpass(&mut self, f0: T, bw: T, slope: T, fs: T) {
        self.set_params(f0, T::zero(), bw, slope, true, fs);
        if self.odd_order {
            self.coeffs.iir1 = IIR1Coefficients::highpass(f0, fs);
            if self.u_slope <= 1 {
                self.update_coeffs(self.coeffs);
                self.process = FilterBand::process_iir1_only;
                return;
            }
            self.process = FilterBand::process_odd_order_cascade;
        } else {
            self.process = FilterBand::process_even_order_cascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        for i in self.start_pole..self.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
            self.coeffs.iir2[i] = IIR2Coefficients::highpass(f0, q_value * q_offset, fs);
        }
        self.update_coeffs(self.coeffs);
    }

    pub fn lowshelf(&mut self, f0: T, gain: T, bw: T, slope: T, fs: T) {
        self.set_params(f0, gain, bw, slope, true, fs);
        let mut partial_gain = gain / self.u_slope.into();
        if self.odd_order {
            self.coeffs.iir1 = IIR1Coefficients::lowshelf(f0, partial_gain, fs);
            if self.u_slope <= 1 {
                self.update_coeffs(self.coeffs);
                self.process = FilterBand::process_iir1_only;
                return;
            }
            self.process = FilterBand::process_odd_order_cascade;
        } else {
            self.process = FilterBand::process_even_order_cascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        partial_gain = partial_gain * Into::<T>::into(2.0);
        for i in self.start_pole..self.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
            self.coeffs.iir2[i] =
                IIR2Coefficients::lowshelf(f0, q_value * q_offset, partial_gain, fs);
        }
        self.update_coeffs(self.coeffs);
    }

    pub fn highshelf(&mut self, f0: T, gain: T, bw: T, slope: T, fs: T) {
        self.set_params(f0, gain, bw, slope, true, fs);
        let mut partial_gain = gain / self.u_slope.into();
        if self.odd_order {
            self.coeffs.iir1 = IIR1Coefficients::highshelf(f0, partial_gain, fs);
            if self.u_slope <= 1 {
                self.update_coeffs(self.coeffs);
                self.process = FilterBand::process_iir1_only;
                return;
            }
            self.process = FilterBand::process_odd_order_cascade;
        } else {
            self.process = FilterBand::process_even_order_cascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        partial_gain = partial_gain * Into::<T>::into(2.0);
        for i in self.start_pole..self.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
            self.coeffs.iir2[i] =
                IIR2Coefficients::highshelf(f0, q_value * q_offset, partial_gain, fs);
        }
        self.update_coeffs(self.coeffs);
    }

    pub fn notch(&mut self, f0: T, gain: T, bw: T, fs: T) {
        let q_value = bw.bw_to_q(f0, fs);
        self.set_params(f0, gain, bw, Into::<T>::into(2.0), false, fs);
        self.coeffs.iir2[0] = IIR2Coefficients::notch(f0, q_value, fs);
        self.process = FilterBand::process_iir2_only;
        self.update_coeffs(self.coeffs);
    }

    pub fn bandpass(&mut self, f0: T, gain: T, bw: T, fs: T) {
        let q_value = bw.bw_to_q(f0, fs);
        self.set_params(f0, gain, bw, Into::<T>::into(2.0), false, fs);
        self.coeffs.iir2[0] = IIR2Coefficients::bandpass(f0, q_value, fs);
        self.process = FilterBand::process_iir2_only;
        self.update_coeffs(self.coeffs);
    }

    pub fn allpass(&mut self, f0: T, bw: T, slope: T, fs: T) {
        self.set_params(f0, T::zero(), bw, slope, true, fs);
        if self.odd_order {
            self.coeffs.iir1 = IIR1Coefficients::allpass(f0, fs);
            if self.u_slope <= 1 {
                self.update_coeffs(self.coeffs);
                self.process = FilterBand::process_iir1_only;
                return;
            }
            self.process = FilterBand::process_odd_order_cascade;
        } else {
            self.process = FilterBand::process_even_order_cascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        for i in self.start_pole..self.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
            self.coeffs.iir2[i] = IIR2Coefficients::allpass(f0, q_value * q_offset, fs);
        }
        self.update_coeffs(self.coeffs);
    }

    pub fn bell(&mut self, f0: T, gain: T, bw: T, fs: T) {
        let q_value = bw.bw_to_q(f0, fs);
        self.set_params(f0, gain, bw, Into::<T>::into(2.0), false, fs);
        self.coeffs.iir2[0] = IIR2Coefficients::bell(f0, q_value, gain, fs);
        self.process = FilterBand::process_iir2_only;
        self.update_coeffs(self.coeffs);
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
        let mut right: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();

        let sample_rate = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let bandwidth = 1.0;
        let slope = 4.0;

        let mut filter_left = FilterBand::new(sample_rate);
        let mut filter_right = FilterBand::new(sample_rate);

        filter_left.highshelf(f0, gain, bandwidth, slope, sample_rate);
        filter_right.mimic_band(&filter_left);

        for i in 0..1000 {
            left[i] = (filter_left.process)(&mut filter_left, left[i]);
            right[i] = (filter_right.process)(&mut filter_right, right[i]);
        }
    }
}
