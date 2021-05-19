use crate::{
    filter_band::{FilterBandCoefficients, ProcessType},
    first_order_iir_wide::{WideF64IIR1, WideF64IIR1Coefficients},
    second_order_iir_wide::{WideF64IIR2, WideF64IIR2Coefficients},
    units::FP,
    MAX_POLE_COUNT,
};
use wide::f64x4;
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
