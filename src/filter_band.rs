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
    units::{Units, ZSample},
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
            let mut y = self.iir1.get_bode_sample(z.pow1);
            assert!(self.iir2.len() >= self.iir2_cascade_count);
            for i in 0..self.iir2_cascade_count {
                y = y * self.iir2[i].get_bode_sample(z);
            }
            y
        } else {
            let mut y = self.iir2[0].get_bode_sample(z);
            assert!(self.iir2.len() >= self.iir2_cascade_count);
            for i in 1..self.iir2_cascade_count {
                y = y * self.iir2[i].get_bode_sample(z);
            }
            y
        }
    }

    pub fn lowpass(
        cutoff_hz: T,
        bandwidth_oct: T,
        order: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            cutoff_hz,
            bandwidth_oct,
            order,
            T::N0,
            sample_rate_hz,
            IIR1Coefficients::lowpass,
            IIR2Coefficients::lowpass,
        )
    }

    pub fn highpass(
        cutoff_hz: T,
        bandwidth_oct: T,
        order: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            cutoff_hz,
            bandwidth_oct,
            order,
            T::N0,
            sample_rate_hz,
            IIR1Coefficients::highpass,
            IIR2Coefficients::highpass,
        )
    }

    pub fn allpass(
        cutoff_hz: T,
        bandwidth_oct: T,
        order: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            cutoff_hz,
            bandwidth_oct,
            order,
            T::N0,
            sample_rate_hz,
            IIR1Coefficients::allpass,
            IIR2Coefficients::allpass,
        )
    }

    pub fn lowshelf(
        cutoff_hz: T,
        gain_db: T,
        bandwidth_oct: T,
        order: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            cutoff_hz,
            bandwidth_oct,
            order,
            gain_db,
            sample_rate_hz,
            IIR1Coefficients::lowshelf,
            IIR2Coefficients::lowshelf,
        )
    }

    pub fn highshelf(
        cutoff_hz: T,
        gain_db: T,
        bandwidth_oct: T,
        order: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        FilterBandCoefficients::filter_type_1(
            cutoff_hz,
            bandwidth_oct,
            order,
            gain_db,
            sample_rate_hz,
            IIR1Coefficients::highshelf,
            IIR2Coefficients::highshelf,
        )
    }

    pub fn filter_type_1(
        cutoff_hz: T,
        bandwidth_oct: T,
        order: T,
        gain_db: T,
        sample_rate_hz: T,
        iir1_coeff_func: fn(T, T, T) -> IIR1Coefficients<T>,
        iir2_coeff_func: fn(T, T, T, T) -> IIR2Coefficients<T>,
    ) -> FilterBandCoefficients<T> {
        let order = order.floor();
        let odd_order = order % T::N2;
        let iir1_enabled = odd_order == T::N1;
        let mut partial_gain = gain_db / order;
        let mut iir1 = IIR1Coefficients::empty();
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        let mut process = ProcessType::ProcessIIR1Only;
        if iir1_enabled {
            iir1 = (iir1_coeff_func)(cutoff_hz, partial_gain, sample_rate_hz);
            if order <= T::N1 {
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
        partial_gain = partial_gain * T::N2;
        let q_value = bandwidth_oct.bandwidth_to_q(cutoff_hz, sample_rate_hz);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q
        let iir2_cascade_count = NumCast::from((order - odd_order) / T::N2).unwrap();
        let odd_order_usize: usize = NumCast::from(odd_order).unwrap();
        let order_usize: usize = NumCast::from(order).unwrap();
        assert!(order_usize < T::BUTTERWORTH.len());
        assert!(odd_order_usize < T::BUTTERWORTH[0].len());
        for i in 0usize..iir2_cascade_count {
            //let q_value: T = butterworth_cascade_q(order_usize, i + odd_order_usize);
            let q_value = T::BUTTERWORTH[order_usize][i + odd_order_usize];
            iir2[i] =
                (iir2_coeff_func)(cutoff_hz, partial_gain, q_value * q_offset, sample_rate_hz);
        }
        FilterBandCoefficients {
            iir1,
            iir2,
            process,
            iir2_cascade_count,
            iir1_enabled,
        }
    }

    pub fn notch(
        cutoff_hz: T,
        _gain_db: T,
        bandwidth_oct: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        iir2[0] = IIR2Coefficients::notch(
            cutoff_hz,
            T::N0,
            bandwidth_oct.bandwidth_to_q(cutoff_hz, sample_rate_hz),
            sample_rate_hz,
        );
        FilterBandCoefficients {
            iir1: IIR1Coefficients::empty(),
            iir2,
            process: ProcessType::ProcessIIR2Only,
            iir2_cascade_count: 0,
            iir1_enabled: false,
        }
    }

    pub fn bandpass(
        cutoff_hz: T,
        _gain_db: T,
        bandwidth_oct: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        iir2[0] = IIR2Coefficients::bandpass(
            cutoff_hz,
            T::N0,
            bandwidth_oct.bandwidth_to_q(cutoff_hz, sample_rate_hz),
            sample_rate_hz,
        );
        FilterBandCoefficients {
            iir1: IIR1Coefficients::empty(),
            iir2,
            process: ProcessType::ProcessIIR2Only,
            iir2_cascade_count: 0,
            iir1_enabled: false,
        }
    }

    pub fn bell(
        cutoff_hz: T,
        gain_db: T,
        bandwidth_oct: T,
        sample_rate_hz: T,
    ) -> FilterBandCoefficients<T> {
        let mut iir2 = [IIR2Coefficients::empty(); MAX_POLE_COUNT];
        iir2[0] = IIR2Coefficients::bell(
            cutoff_hz,
            gain_db,
            bandwidth_oct.bandwidth_to_q(cutoff_hz, sample_rate_hz),
            sample_rate_hz,
        );
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

    pub fn process_iir1_only(&mut self, input_sample: T) -> T {
        self.iir1.process(input_sample)
    }

    pub fn process_iir2_only(&mut self, input_sample: T) -> T {
        self.iir2[0].process(input_sample)
    }

    pub fn process_even_order_cascade(&mut self, input_sample: T) -> T {
        assert!(self.iir2.len() >= self.iir2_cascade_count);
        let mut input_sample = input_sample;
        for i in 0..self.iir2_cascade_count {
            input_sample = self.iir2[i].process(input_sample);
        }
        input_sample
    }

    pub fn process_odd_order_cascade(&mut self, input_sample: T) -> T {
        assert!(self.iir2.len() >= self.iir2_cascade_count);
        let mut input_sample = self.iir1.process(input_sample);
        for i in 0..self.iir2_cascade_count {
            input_sample = self.iir2[i].process(input_sample);
        }
        input_sample
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
    fn test_filter_band() {
        let mut left: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();
        let mut right: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();

        let fs = 48000.0;
        let cutoff_hz = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let order = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(cutoff_hz, gain, width, order, fs);

        let mut filter_left = FilterBand::from(&coeffs);
        let mut filter_right = FilterBand::from(&coeffs);

        for i in 0..1000 {
            left[i] = (filter_left.process)(&mut filter_left, left[i]);
            right[i] = (filter_right.process)(&mut filter_right, right[i]);
        }

        dbg!(left[500], right[500]);
    }
}
