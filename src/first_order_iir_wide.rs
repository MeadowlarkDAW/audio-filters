use crate::{first_order_iir::IIR1Coefficients, units::FP, wide_units::WIDE};

#[derive(Copy, Clone, Debug)]
pub struct WideIIR1Coefficients<T: WIDE> {
    pub a: T,
    pub g: T,
    pub a1: T,
    pub m0: T,
    pub m1: T,
    pub fs: T,
}

impl<T: WIDE> WideIIR1Coefficients<T> {
    pub fn from<A: FP>(coeffs: IIR1Coefficients<A>) -> WideIIR1Coefficients<T> {
        let a = T::from_w(coeffs.a);
        let g = T::from_w(coeffs.g);
        let a1 = T::from_w(coeffs.a1);
        let m0 = T::from_w(coeffs.m0);
        let m1 = T::from_w(coeffs.m1);
        let fs = T::from_w(coeffs.fs);
        WideIIR1Coefficients {
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
pub struct WideIIR1<T: WIDE> {
    ic1eq: T,
    pub coeffs: WideIIR1Coefficients<T>,
}

impl<T: WIDE> WideIIR1<T> {
    pub fn new(coefficients: WideIIR1Coefficients<T>) -> Self {
        WideIIR1 {
            ic1eq: T::ZERO,
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: T) -> T {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input + self.coeffs.m1 * v2
    }

    pub fn process_partial(&mut self, input: T) -> T {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: WideIIR1Coefficients<T>) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
mod tests {
    use wide::f64x4;

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
        let coeffs = WideIIR1Coefficients::from(coeffs);

        let mut filter_left = WideIIR1::new(coeffs);

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
