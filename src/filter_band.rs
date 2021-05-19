use core::u8;

use crate::{
    first_order_iir::{WideF64IIR1, WideF64IIR1Coefficients},
    second_order_iir::{WideF64IIR2, WideF64IIR2Coefficients},
    units::FP,
};
use num_complex::Complex;

use num_traits::NumCast;

use wide::f64x4;

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
    pub process: ProcessType,
    pub iir2_cascade_count: usize,
    pub iir1_enabled: bool,
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

    pub fn lowpass(f0: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            f0,
            bw,
            slope,
            T::zero(),
            fs,
            IIR1Coefficients::lowpass,
            IIR2Coefficients::lowpass,
        )
    }

    pub fn highpass(f0: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            f0,
            bw,
            slope,
            T::zero(),
            fs,
            IIR1Coefficients::highpass,
            IIR2Coefficients::highpass,
        )
    }

    pub fn allpass(f0: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            f0,
            bw,
            slope,
            T::zero(),
            fs,
            IIR1Coefficients::allpass,
            IIR2Coefficients::allpass,
        )
    }

    pub fn lowshelf(f0: T, gain: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            f0,
            bw,
            slope,
            gain,
            fs,
            IIR1Coefficients::lowshelf,
            IIR2Coefficients::lowshelf,
        )
    }

    pub fn highshelf(f0: T, gain: T, bw: T, slope: T, fs: T) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            f0,
            bw,
            slope,
            gain,
            fs,
            IIR1Coefficients::highshelf,
            IIR2Coefficients::highshelf,
        )
    }

    pub fn filter_type_1(
        f0: T,
        bw: T,
        slope: T,
        gain: T,
        fs: T,
        iir1_func: fn(T, T, T) -> IIR1Coefficients<T>,
        iir2_func: fn(T, T, T, T) -> IIR2Coefficients<T>,
    ) -> FilterBandCoefficients<T> {
        let u_slope: usize = NumCast::from(slope).unwrap();
        let odd_order = u_slope & 1usize == 1usize;
        let iir1_enabled = odd_order;
        let start_pole = iir1_enabled as usize;
        let mut partial_gain = gain / NumCast::from(u_slope).unwrap();
        let mut iir1 = IIR1Coefficients::empty();
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        let mut process = ProcessType::ProcessIIR1Only;
        if odd_order {
            iir1 = (iir1_func)(f0, partial_gain, fs);
            if u_slope <= 1 {
                return FilterBandCoefficients {
                    iir1,
                    iir2,
                    process,
                    iir2_cascade_count: 0,
                    iir1_enabled,
                };
            }
            process = ProcessType::ProcessOddOrderCascade;
        } else {
            process = ProcessType::ProcessEvenOrderCascade;
        }
        partial_gain = partial_gain * Into::<T>::into(2.0);
        let q_value = bw.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        let iir2_cascade_count = ((u_slope as usize - start_pole) / 2usize) as usize;
        for i in 0..iir2_cascade_count {
            let q_value: T = butterworth_cascade_q(u_slope as u8, i as u8 + start_pole as u8);
            iir2[i] = (iir2_func)(f0, q_value * q_offset, partial_gain, fs);
        }
        FilterBandCoefficients {
            iir1,
            iir2,
            process,
            iir2_cascade_count,
            iir1_enabled,
        }
    }

    pub fn notch(f0: T, _gain: T, bw: T, fs: T) -> FilterBandCoefficients<T> {
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        iir2[0] = IIR2Coefficients::notch(f0, bw.bw_to_q(f0, fs), T::zero(), fs);
        FilterBandCoefficients {
            iir1: IIR1Coefficients::empty(),
            iir2,
            process: ProcessType::ProcessIIR2Only,
            iir2_cascade_count: 0,
            iir1_enabled: false,
        }
    }

    pub fn bandpass(f0: T, _gain: T, bw: T, fs: T) -> FilterBandCoefficients<T> {
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        iir2[0] = IIR2Coefficients::bandpass(f0, bw.bw_to_q(f0, fs), T::zero(), fs);
        FilterBandCoefficients {
            iir1: IIR1Coefficients::empty(),
            iir2,
            process: ProcessType::ProcessIIR2Only,
            iir2_cascade_count: 0,
            iir1_enabled: false,
        }
    }

    pub fn bell(f0: T, gain: T, bw: T, fs: T) -> FilterBandCoefficients<T> {
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        iir2[0] = IIR2Coefficients::bell(f0, bw.bw_to_q(f0, fs), gain, fs);
        FilterBandCoefficients {
            iir1: IIR1Coefficients::empty(),
            iir2,
            process: ProcessType::ProcessIIR2Only,
            iir2_cascade_count: 0,
            iir1_enabled: false,
        }
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

#[derive(Copy, Clone, Debug)]
pub struct WideF64FilterBandCoefficients {
    pub iir1: WideF64IIR1Coefficients,
    pub iir2: [WideF64IIR2Coefficients; MAX_POLE_COUNT],
    pub process: ProcessType,
    pub iir2_cascade_count: usize,
    pub iir1_enabled: bool,
}

impl WideF64FilterBandCoefficients {
    pub fn from<T: FP>(coeffs: FilterBandCoefficients<T>) -> WideF64FilterBandCoefficients {
        let mut iir2_cascade = WideF64IIR2Coefficients::empty_cascade();
        for (iir2, in_iir2) in iir2_cascade.iter_mut().zip(&coeffs.iir2) {
            *iir2 = WideF64IIR2Coefficients::from(*in_iir2);
        }
        WideF64FilterBandCoefficients {
            iir1: WideF64IIR1Coefficients::from(coeffs.iir1),
            iir2: iir2_cascade,
            process: coeffs.process,
            iir2_cascade_count: coeffs.iir2_cascade_count,
            iir1_enabled: coeffs.iir1_enabled,
        }
    }
}

#[derive(Copy, Clone)]
pub struct WideF64FilterBand {
    iir1: WideF64IIR1,
    iir2: [WideF64IIR2; MAX_POLE_COUNT],
    iir2_cascade_count: usize,
    pub process: fn(&mut Self, f64x4) -> f64x4,
}

impl WideF64FilterBand {
    pub fn from(coeffs: &WideF64FilterBandCoefficients) -> WideF64FilterBand {
        WideF64FilterBand {
            iir1: WideF64IIR1::new(coeffs.iir1),
            iir2: [WideF64IIR2::new(coeffs.iir2[0]); MAX_POLE_COUNT],
            iir2_cascade_count: coeffs.iir2_cascade_count,
            process: WideF64FilterBand::get_process(coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, x: f64x4) -> f64x4 {
        self.iir1.process(x)
    }

    pub fn process_iir2_only(&mut self, x: f64x4) -> f64x4 {
        self.iir2[0].process(x)
    }

    pub fn process_even_order_cascade(&mut self, x: f64x4) -> f64x4 {
        let mut x = x;
        for i in 0..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn process_odd_order_cascade(&mut self, x: f64x4) -> f64x4 {
        let mut x = self.iir1.process(x);
        for i in 0..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, f64x4) -> f64x4 {
        match process_type {
            ProcessType::ProcessIIR1Only => WideF64FilterBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => WideF64FilterBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => WideF64FilterBand::process_even_order_cascade,
            ProcessType::ProcessOddOrderCascade => WideF64FilterBand::process_odd_order_cascade,
        }
    }

    pub fn update(&mut self, coeffs: &WideF64FilterBandCoefficients) {
        for (filter, coeff) in self.iir2.iter_mut().zip(coeffs.iir2.iter()) {
            filter.update_coefficients(*coeff)
        }
        self.iir1.update_coefficients(coeffs.iir1);
        self.iir2_cascade_count = coeffs.iir2_cascade_count;
        self.process = WideF64FilterBand::get_process(coeffs.process);
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

        let mut filter_left = FilterBand::from(&coeffs);
        let mut filter_right = FilterBand::from(&coeffs);

        for i in 0..1000 {
            left[i] = (filter_left.process)(&mut filter_left, left[i]);
            right[i] = (filter_right.process)(&mut filter_right, right[i]);
        }

        dbg!(left[500], right[500]);
    }

    fn rand64(x: f64) -> f64 {
        ((x * 12.98983123).sin() * 43758.545345345).fract()
    }

    #[test]
    fn test_wide_filter_band() {
        let mut ch1: Vec<f64> = (0..1000).map(|x| rand64(x as f64)).collect();
        let mut ch2: Vec<f64> = (1000..2000).map(|x| rand64(x as f64)).collect();
        let mut ch3: Vec<f64> = (2000..3000).map(|x| rand64(x as f64)).collect();
        let mut ch4: Vec<f64> = (3000..4000).map(|x| rand64(x as f64)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let slope = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);
        let coeffs = WideF64FilterBandCoefficients::from(coeffs);

        let mut filter = WideF64FilterBand::from(&coeffs);

        for i in 0..1000 {
            let output: [f64; 4] =
                (filter.process)(&mut filter, f64x4::from([ch1[i], ch2[i], ch3[i], ch4[i]])).into();
            ch1[i] = output[0];
            ch2[i] = output[1];
            ch3[i] = output[2];
            ch4[i] = output[3];
        }

        println!("{} {} {} {}", ch1[500], ch2[500], ch3[500], ch4[500])
    }
}
