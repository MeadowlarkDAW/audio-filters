use crate::units::FP;
use num_complex::Complex;

use num_traits::NumCast;

#[derive(Clone, Copy, Debug)]
pub enum ProcessType {
    ProcessIIR1Only,
    ProcessIIR2Only,
    ProcessEvenOrderCascade,
    ProcessOddOrderCascade,
}

use crate::{
    first_order_iir::{IIR1Coefficients, IIR1},
    second_order_iir::{IIR2Coefficients, IIR2},
    units::{butterworth_cascade_q, Units, ZSample},
    MAX_POLE_COUNT,
};

#[derive(Copy, Clone, Debug)]
pub struct FilterBandCoefficients<T: FP> {
    pub iir1: IIR1Coefficients<T>,
    pub iir2: [IIR2Coefficients<T>; MAX_POLE_COUNT],
    pub iir1_enabled: bool,
    pub f0: T,
    pub gain: T,
    pub bw_value: T,
    pub slope: T,
    pub u_slope: u8,
    pub odd_order: bool,
    pub start_pole: usize,
    pub iir2_cascade_count: usize,
    pub sample_rate: T,
    pub process: ProcessType,
}

impl<T: FP> FilterBandCoefficients<T> {
    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        if self.iir1_enabled {
            let mut y = self.iir1.get_bode_sample(z.z);
            for i in 0..self.iir2_cascade_count {
                y = y * self.iir2[i].get_bode_sample(z);
            }
            y
        } else {
            let mut y = self.iir2[0].get_bode_sample(z);
            for i in 1..self.iir2_cascade_count {
                y = y * self.iir2[i].get_bode_sample(z);
            }
            y
        }
    }

    fn init(
        f0: T,
        gain: T,
        slope: T,
        bw: T,
        iir1_possible: bool,
        fs: T,
    ) -> FilterBandCoefficients<T> {
        let iir2_coeffs = IIR2Coefficients::bell(T::zero(), T::one(), T::zero(), fs);
        let iir1_coeffs = IIR1Coefficients::<T>::lowpass(1000.0.into(), fs);

        let u_slope = NumCast::from(slope).unwrap();

        let odd_order = u_slope & 1 == 1;
        let iir1_enabled = iir1_possible && odd_order;
        let start_pole = iir1_enabled as usize;

        FilterBandCoefficients {
            iir1: iir1_coeffs,
            iir2: [iir2_coeffs; MAX_POLE_COUNT],
            iir1_enabled,
            f0,
            gain,
            bw_value: bw,
            slope,
            odd_order,
            u_slope,
            start_pole,
            iir2_cascade_count: ((u_slope as usize - start_pole) / 2) as usize,
            sample_rate: fs,
            process: ProcessType::ProcessIIR1Only,
        }
    }

    pub fn lowpass(f0: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        let mut co = FilterBandCoefficients::init(f0, T::zero(), slope, bw, true, fs);
        if co.odd_order {
            co.iir1 = IIR1Coefficients::lowpass(f0, fs);
            if co.u_slope <= 1 {
                co.process = ProcessType::ProcessIIR1Only;
                return co;
            }
            co.process = ProcessType::ProcessOddOrderCascade;
        } else {
            co.process = ProcessType::ProcessEvenOrderCascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        for i in 0..co.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(co.u_slope, i as u8 + co.start_pole as u8);
            co.iir2[i] = IIR2Coefficients::lowpass(f0, q_value * q_offset, fs);
        }
        co
    }

    pub fn highpass(f0: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        let mut co = FilterBandCoefficients::init(f0, T::zero(), slope, bw, true, fs);
        if co.odd_order {
            co.iir1 = IIR1Coefficients::highpass(f0, fs);
            if co.u_slope <= 1 {
                co.process = ProcessType::ProcessIIR1Only;
                return co;
            }
            co.process = ProcessType::ProcessOddOrderCascade;
        } else {
            co.process = ProcessType::ProcessEvenOrderCascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        for i in 0..co.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(co.u_slope, i as u8 + co.start_pole as u8);
            co.iir2[i] = IIR2Coefficients::highpass(f0, q_value * q_offset, fs);
        }
        co
    }

    pub fn lowshelf(f0: T, gain: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        let mut co = FilterBandCoefficients::init(f0, gain, slope, bw, true, fs);
        let mut partial_gain = gain / co.u_slope.into();
        if co.odd_order {
            co.iir1 = IIR1Coefficients::lowshelf(f0, partial_gain, fs);
            if co.u_slope <= 1 {
                co.process = ProcessType::ProcessIIR1Only;
                return co;
            }
            co.process = ProcessType::ProcessOddOrderCascade;
        } else {
            co.process = ProcessType::ProcessEvenOrderCascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        partial_gain = partial_gain * Into::<T>::into(2.0);
        for i in 0..co.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(co.u_slope, i as u8 + co.start_pole as u8);
            co.iir2[i] = IIR2Coefficients::lowshelf(f0, q_value * q_offset, partial_gain, fs);
        }
        co
    }

    pub fn highshelf(f0: T, gain: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        let mut co = FilterBandCoefficients::init(f0, gain, slope, bw, true, fs);
        let mut partial_gain = gain / co.u_slope.into();
        if co.odd_order {
            co.iir1 = IIR1Coefficients::highshelf(f0, partial_gain, fs);
            if co.u_slope <= 1 {
                co.process = ProcessType::ProcessIIR1Only;
                return co;
            }
            co.process = ProcessType::ProcessOddOrderCascade;
        } else {
            co.process = ProcessType::ProcessEvenOrderCascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        partial_gain = partial_gain * Into::<T>::into(2.0);
        for i in 0..co.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(co.u_slope, i as u8 + co.start_pole as u8);
            co.iir2[i] = IIR2Coefficients::highshelf(f0, q_value * q_offset, partial_gain, fs);
        }
        co
    }

    pub fn notch(f0: T, gain: T, bw: T, fs: T) -> FilterBandCoefficients<T> {
        let q_value = bw.bw_to_q(f0, fs);
        let mut co = FilterBandCoefficients::init(f0, gain, Into::<T>::into(2.0), bw, false, fs);
        co.iir2[0] = IIR2Coefficients::notch(f0, q_value, fs);
        co.process = ProcessType::ProcessIIR2Only;
        co
    }

    pub fn bandpass(f0: T, gain: T, bw: T, fs: T) -> FilterBandCoefficients<T> {
        let q_value = bw.bw_to_q(f0, fs);
        let mut co = FilterBandCoefficients::init(f0, gain, Into::<T>::into(2.0), bw, false, fs);
        co.iir2[0] = IIR2Coefficients::bandpass(f0, q_value, fs);
        co.process = ProcessType::ProcessIIR2Only;
        co
    }

    pub fn allpass(f0: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        let mut co = FilterBandCoefficients::init(f0, T::zero(), slope, bw, true, fs);
        if co.odd_order {
            co.iir1 = IIR1Coefficients::allpass(f0, fs);
            if co.u_slope <= 1 {
                co.process = ProcessType::ProcessIIR1Only;
                return co;
            }
            co.process = ProcessType::ProcessOddOrderCascade;
        } else {
            co.process = ProcessType::ProcessEvenOrderCascade;
        }
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        for i in 0..co.iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(co.u_slope, i as u8 + co.start_pole as u8);
            co.iir2[i] = IIR2Coefficients::allpass(f0, q_value * q_offset, fs);
        }
        co
    }

    pub fn bell(f0: T, gain: T, bw: T, fs: T) -> FilterBandCoefficients<T> {
        let q_value = bw.bw_to_q(f0, fs);
        let mut co = FilterBandCoefficients::init(f0, gain, Into::<T>::into(2.0), bw, false, fs);
        co.iir2[0] = IIR2Coefficients::bell(f0, q_value, gain, fs);
        co.process = ProcessType::ProcessIIR2Only;
        co
    }
}

#[derive(Copy, Clone)]
pub struct FilterBand<T: FP> {
    iir1: IIR1<T>,
    iir2: [IIR2<T>; MAX_POLE_COUNT],
    iir2_cascade_count: usize,
    pub process: fn(&mut Self, T) -> T,
}

impl<T: FP> FilterBand<T> {
    pub fn from(coeffs: &FilterBandCoefficients<T>) -> FilterBand<T> {
        FilterBand {
            iir1: IIR1::<T>::new(coeffs.iir1),
            iir2: [IIR2::<T>::new(coeffs.iir2[0]); MAX_POLE_COUNT],
            iir2_cascade_count: coeffs.iir2_cascade_count,
            process: FilterBand::get_process(coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, x: T) -> T {
        self.iir1.process(x)
    }

    pub fn process_iir2_only(&mut self, x: T) -> T {
        self.iir2[0].process(x)
    }

    pub fn process_even_order_cascade(&mut self, x: T) -> T {
        let mut x = x;
        for i in 0..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn process_odd_order_cascade(&mut self, x: T) -> T {
        let mut x = self.iir1.process(x);
        for i in 0..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, T) -> T {
        match process_type {
            ProcessType::ProcessIIR1Only => FilterBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => FilterBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => FilterBand::process_even_order_cascade,
            ProcessType::ProcessOddOrderCascade => FilterBand::process_odd_order_cascade,
        }
    }

    pub fn update(&mut self, coeffs: &FilterBandCoefficients<T>) {
        for (filter, coeff) in self.iir2.iter_mut().zip(coeffs.iir2.iter()) {
            filter.update_coefficients(*coeff)
        }
        self.iir1.update_coefficients(coeffs.iir1);
        self.iir2_cascade_count = coeffs.iir2_cascade_count;
        self.process = FilterBand::get_process(coeffs.process);
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

        let fs = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let slope = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);

        let mut filter_left = FilterBand::from(&coeffs);
        let mut filter_right = FilterBand::from(&coeffs);

        for i in 0..1000 {
            left[i] = (filter_left.process)(&mut filter_left, left[i]);
            right[i] = (filter_right.process)(&mut filter_right, right[i]);
        }

        dbg!(left[500], right[500]);
    }
}
