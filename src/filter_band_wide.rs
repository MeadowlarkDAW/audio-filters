use crate::{
    filter_band::{FilterBandCoefficients, ProcessType},
    first_order_iir_wide::{WideIIR1, WideIIR1Coefficients},
    second_order_iir_wide::{WideIIR2, WideIIR2Coefficients},
    units::FP,
    wide_units::WIDE,
    MAX_POLE_COUNT,
};

#[derive(Copy, Clone, Debug)]
pub struct WideFilterBandCoefficients<T: WIDE> {
    pub iir1: WideIIR1Coefficients<T>,
    pub iir2: [WideIIR2Coefficients<T>; MAX_POLE_COUNT],
    pub process: ProcessType,
    pub iir2_cascade_count: usize,
    pub iir1_enabled: bool,
}

impl<T: WIDE> WideFilterBandCoefficients<T> {
    pub fn from<A: FP>(coeffs: FilterBandCoefficients<A>) -> WideFilterBandCoefficients<T> {
        let mut iir2_cascade = WideIIR2Coefficients::empty_cascade();
        for (iir2, in_iir2) in iir2_cascade.iter_mut().zip(&coeffs.iir2) {
            *iir2 = WideIIR2Coefficients::from(*in_iir2);
        }
        WideFilterBandCoefficients {
            iir1: WideIIR1Coefficients::from(coeffs.iir1),
            iir2: iir2_cascade,
            process: coeffs.process,
            iir2_cascade_count: coeffs.iir2_cascade_count,
            iir1_enabled: coeffs.iir1_enabled,
        }
    }
}

#[derive(Copy, Clone)]
pub struct WideFilterBand<T: WIDE> {
    iir1: WideIIR1<T>,
    iir2: [WideIIR2<T>; MAX_POLE_COUNT],
    iir2_cascade_count: usize,
    pub process: fn(&mut Self, T) -> T,
}

impl<T: WIDE> WideFilterBand<T> {
    pub fn from(coeffs: &WideFilterBandCoefficients<T>) -> WideFilterBand<T> {
        WideFilterBand {
            iir1: WideIIR1::new(coeffs.iir1),
            iir2: [WideIIR2::new(coeffs.iir2[0]); MAX_POLE_COUNT],
            iir2_cascade_count: coeffs.iir2_cascade_count,
            process: WideFilterBand::get_process(coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, x: T) -> T {
        self.iir1.process(x)
    }

    pub fn process_iir2_only(&mut self, x: T) -> T {
        self.iir2[0].process(x)
    }

    pub fn process_even_order_cascade(&mut self, x: T) -> T {
        assert!(self.iir2.len() >= self.iir2_cascade_count);
        let mut x = x;
        for i in 0..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn process_odd_order_cascade(&mut self, x: T) -> T {
        assert!(self.iir2.len() >= self.iir2_cascade_count);
        let mut x = self.iir1.process(x);
        for i in 0..self.iir2_cascade_count {
            x = self.iir2[i].process(x);
        }
        x
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, T) -> T {
        match process_type {
            ProcessType::ProcessIIR1Only => WideFilterBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => WideFilterBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => WideFilterBand::process_even_order_cascade,
            ProcessType::ProcessOddOrderCascade => WideFilterBand::process_odd_order_cascade,
        }
    }

    pub fn update(&mut self, coeffs: &WideFilterBandCoefficients<T>) {
        for (filter, coeff) in self.iir2.iter_mut().zip(coeffs.iir2.iter()) {
            filter.update_coefficients(*coeff)
        }
        self.iir1.update_coefficients(coeffs.iir1);
        self.iir2_cascade_count = coeffs.iir2_cascade_count;
        self.process = WideFilterBand::get_process(coeffs.process);
    }
}

#[cfg(test)]
mod tests {
    use wide::f32x4;
    use wide::f32x8;
    use wide::f64x2;
    use wide::f64x4;

    use super::*;

    fn rand64(x: f64) -> f64 {
        ((x * 12.989846024374758).sin() * 43758.545347294991945).fract()
    }

    fn rand32(x: f32) -> f32 {
        ((x * 12.989846024374758).sin() * 43758.545347294991945).fract()
    }

    #[test]
    fn test_widef64x4_filter_band() {
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
        let coeffs = WideFilterBandCoefficients::from(coeffs);

        let mut filter = WideFilterBand::from(&coeffs);

        for i in 0..1000 {
            let output: [f64; 4] =
                (filter.process)(&mut filter, f64x4::from([ch1[i], ch2[i], ch3[i], ch4[i]])).into();
            ch1[i] = output[0];
            ch2[i] = output[1];
            ch3[i] = output[2];
            ch4[i] = output[3];
        }

        assert_eq!(
            [
                -1.8667512705219382,
                2.0001785972507573,
                1.2690387167880355,
                -1.4693735156714782
            ],
            [ch1[500], ch2[500], ch3[500], ch4[500]]
        );
    }

    #[test]
    fn test_widef64x2_filter_band() {
        let mut ch1: Vec<f64> = (0..1000).map(|x| rand64(x as f64)).collect();
        let mut ch2: Vec<f64> = (1000..2000).map(|x| rand64(x as f64)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let slope = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);
        let coeffs = WideFilterBandCoefficients::from(coeffs);

        let mut filter = WideFilterBand::from(&coeffs);

        for i in 0..1000 {
            let output: [f64; 2] =
                (filter.process)(&mut filter, f64x2::from([ch1[i], ch2[i]])).into();
            ch1[i] = output[0];
            ch2[i] = output[1];
        }

        assert_eq!(
            [-1.8667512705219382, 2.0001785972507573],
            [ch1[500], ch2[500]]
        );
    }

    #[test]
    fn test_widef32x8_filter_band() {
        let mut ch1: Vec<f32> = (0..1000).map(|x| rand32(x as f32)).collect();
        let mut ch2: Vec<f32> = (1000..2000).map(|x| rand32(x as f32)).collect();
        let mut ch3: Vec<f32> = (2000..3000).map(|x| rand32(x as f32)).collect();
        let mut ch4: Vec<f32> = (3000..4000).map(|x| rand32(x as f32)).collect();
        let mut ch5: Vec<f32> = (4000..5000).map(|x| rand32(x as f32)).collect();
        let mut ch6: Vec<f32> = (5000..6000).map(|x| rand32(x as f32)).collect();
        let mut ch7: Vec<f32> = (6000..7000).map(|x| rand32(x as f32)).collect();
        let mut ch8: Vec<f32> = (7000..8000).map(|x| rand32(x as f32)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let slope = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);
        let coeffs = WideFilterBandCoefficients::from(coeffs);

        let mut filter = WideFilterBand::from(&coeffs);

        for i in 0..1000 {
            let output: [f32; 8] = (filter.process)(
                &mut filter,
                f32x8::from([
                    ch1[i], ch2[i], ch3[i], ch4[i], ch5[i], ch6[i], ch7[i], ch8[i],
                ]),
            )
            .into();
            ch1[i] = output[0];
            ch2[i] = output[1];
            ch3[i] = output[2];
            ch4[i] = output[3];
            ch5[i] = output[4];
            ch6[i] = output[5];
            ch7[i] = output[6];
            ch8[i] = output[7];
        }

        assert_eq!(
            [
                -0.84111476,
                1.4498048,
                0.048668005,
                -1.60332,
                0.6705888,
                -0.72239006,
                1.1677216,
                1.436455
            ],
            [ch1[500], ch2[500], ch3[500], ch4[500], ch5[500], ch6[500], ch7[500], ch8[500]]
        );
    }

    #[test]
    fn test_widef32x4_filter_band() {
        let mut ch1: Vec<f32> = (0..1000).map(|x| rand32(x as f32)).collect();
        let mut ch2: Vec<f32> = (1000..2000).map(|x| rand32(x as f32)).collect();
        let mut ch3: Vec<f32> = (2000..3000).map(|x| rand32(x as f32)).collect();
        let mut ch4: Vec<f32> = (3000..4000).map(|x| rand32(x as f32)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let gain = 6.0;
        let width = 1.0;
        let slope = 4.0;

        let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);
        let coeffs = WideFilterBandCoefficients::from(coeffs);

        let mut filter = WideFilterBand::from(&coeffs);

        for i in 0..1000 {
            let output: [f32; 4] =
                (filter.process)(&mut filter, f32x4::from([ch1[i], ch2[i], ch3[i], ch4[i]])).into();
            ch1[i] = output[0];
            ch2[i] = output[1];
            ch3[i] = output[2];
            ch4[i] = output[3];
        }

        assert_eq!(
            [-0.84111476, 1.4498048, 0.048668005, -1.60332],
            [ch1[500], ch2[500], ch3[500], ch4[500]]
        );
    }
}
