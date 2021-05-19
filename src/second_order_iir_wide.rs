use crate::{second_order_iir::IIR2Coefficients, units::FP, MAX_POLE_COUNT};
use wide::f64x4;

#[derive(Copy, Clone, Debug)]
pub struct WideF64IIR2Coefficients {
    pub a: f64x4,
    pub g: f64x4,
    pub gpow2: f64x4,
    pub k: f64x4,
    pub a1: f64x4,
    pub a2: f64x4,
    pub a3: f64x4,
    pub m0: f64x4,
    pub m1: f64x4,
    pub m2: f64x4,
    pub fs: f64x4,
}

impl WideF64IIR2Coefficients {
    pub fn from<T: FP>(coeffs: IIR2Coefficients<T>) -> WideF64IIR2Coefficients {
        let a = f64x4::splat(Into::<f64>::into(coeffs.a));
        let g = f64x4::splat(Into::<f64>::into(coeffs.g));
        let k = f64x4::splat(Into::<f64>::into(coeffs.k));
        let a1 = f64x4::splat(Into::<f64>::into(coeffs.a1));
        let a2 = f64x4::splat(Into::<f64>::into(coeffs.a2));
        let a3 = f64x4::splat(Into::<f64>::into(coeffs.a3));
        let m0 = f64x4::splat(Into::<f64>::into(coeffs.m0));
        let m1 = f64x4::splat(Into::<f64>::into(coeffs.m1));
        let m2 = f64x4::splat(Into::<f64>::into(coeffs.m2));
        let fs = f64x4::splat(Into::<f64>::into(coeffs.fs));
        WideF64IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }

    pub const fn empty() -> WideF64IIR2Coefficients {
        WideF64IIR2Coefficients {
            a: f64x4::ZERO,
            g: f64x4::ZERO,
            gpow2: f64x4::ZERO,
            k: f64x4::ZERO,
            a1: f64x4::ZERO,
            a2: f64x4::ZERO,
            a3: f64x4::ZERO,
            m0: f64x4::ZERO,
            m1: f64x4::ZERO,
            m2: f64x4::ZERO,
            fs: f64x4::ZERO,
        }
    }

    pub const fn empty_cascade() -> [WideF64IIR2Coefficients; MAX_POLE_COUNT] {
        [WideF64IIR2Coefficients::empty(); MAX_POLE_COUNT]
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct WideF64IIR2 {
    ic1eq: f64x4,
    ic2eq: f64x4,
    pub coeffs: WideF64IIR2Coefficients,
}

impl WideF64IIR2 {
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: WideF64IIR2Coefficients) -> Self {
        WideF64IIR2 {
            ic1eq: f64x4::ZERO,
            ic2eq: f64x4::ZERO,
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: f64x4) -> f64x4 {
        let two: f64x4 = 2.0.into();
        let v3 = input - self.ic2eq;
        let v1 = self.coeffs.a1 * self.ic1eq + self.coeffs.a2 * v3;
        let v2 = self.ic2eq + self.coeffs.a2 * self.ic1eq + self.coeffs.a3 * v3;
        self.ic1eq = two * v1 - self.ic1eq;
        self.ic2eq = two * v2 - self.ic2eq;

        self.coeffs.m0 * input + self.coeffs.m1 * v1 + self.coeffs.m2 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: WideF64IIR2Coefficients) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
mod tests {
    use crate::units::Units;

    use super::*;

    fn rand(x: f64) -> f64 {
        ((x * 12.98983123).sin() * 43758.545345345).fract()
    }

    #[test]
    fn wide_test() {
        let mut ch1: Vec<f64> = (0..1000).map(|x| rand(x as f64)).collect();
        let mut ch2: Vec<f64> = (1000..2000).map(|x| rand(x as f64)).collect();
        let mut ch3: Vec<f64> = (2000..3000).map(|x| rand(x as f64)).collect();
        let mut ch4: Vec<f64> = (3000..4000).map(|x| rand(x as f64)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let bandwidth = 1.0;

        let coeffs = IIR2Coefficients::lowpass(f0, bandwidth.bw_to_q(f0, fs), 0.0, fs);
        let coeffs = WideF64IIR2Coefficients::from(coeffs);

        let mut filter_left = WideF64IIR2::new(coeffs);

        for i in 0..1000 {
            let output: [f64; 4] = filter_left
                .process(f64x4::from([ch1[i], ch2[i], ch3[i], ch4[i]]))
                .into();
            ch1[i] = output[0];
            ch2[i] = output[1];
            ch3[i] = output[2];
            ch4[i] = output[3];
        }
        println!("{} {} {} {}", ch1[500], ch2[500], ch3[500], ch4[500])
    }
}
