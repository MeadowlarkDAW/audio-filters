use crate::{first_order_iir::IIR1Coefficients, units::FP};
use wide::f64x4;
#[derive(Copy, Clone, Debug)]
pub struct WideF64IIR1Coefficients {
    pub a: f64x4,
    pub g: f64x4,
    pub a1: f64x4,
    pub m0: f64x4,
    pub m1: f64x4,
    pub fs: f64x4,
}

impl WideF64IIR1Coefficients {
    pub fn from<T: FP>(coeffs: IIR1Coefficients<T>) -> WideF64IIR1Coefficients {
        let a = f64x4::splat(Into::<f64>::into(coeffs.a));
        let g = f64x4::splat(Into::<f64>::into(coeffs.g));
        let a1 = f64x4::splat(Into::<f64>::into(coeffs.a1));
        let m0 = f64x4::splat(Into::<f64>::into(coeffs.m0));
        let m1 = f64x4::splat(Into::<f64>::into(coeffs.m1));
        let fs = f64x4::splat(Into::<f64>::into(coeffs.fs));
        WideF64IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct WideF64IIR1 {
    ic1eq: f64x4,
    pub coeffs: WideF64IIR1Coefficients,
}

impl WideF64IIR1 {
    pub fn new(coefficients: WideF64IIR1Coefficients) -> Self {
        WideF64IIR1 {
            ic1eq: f64x4::ZERO,
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: f64x4) -> f64x4 {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input + self.coeffs.m1 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: WideF64IIR1Coefficients) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
mod tests {
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

        let coeffs = IIR1Coefficients::lowpass(f0, 0.0, fs);
        let coeffs = WideF64IIR1Coefficients::from(coeffs);

        let mut filter_left = WideF64IIR1::new(coeffs);

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
