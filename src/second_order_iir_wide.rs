use crate::{second_order_iir::IIR2Coefficients, units::FP, wide_units::WIDE, MAX_POLE_COUNT};

#[derive(Copy, Clone, Debug)]
pub struct WideIIR2Coefficients<T: WIDE> {
    pub a: T,
    pub g: T,
    pub gpow2: T,
    pub k: T,
    pub a1: T,
    pub a2: T,
    pub a3: T,
    pub m0: T,
    pub m1: T,
    pub m2: T,
    pub fs: T,
}

impl<T: WIDE> WideIIR2Coefficients<T> {
    pub fn from<A: FP>(coeffs: IIR2Coefficients<A>) -> WideIIR2Coefficients<T> {
        let a = T::from_w(coeffs.a);
        let g = T::from_w(coeffs.g);
        let k = T::from_w(coeffs.k);
        let a1 = T::from_w(coeffs.a1);
        let a2 = T::from_w(coeffs.a2);
        let a3 = T::from_w(coeffs.a3);
        let m0 = T::from_w(coeffs.m0);
        let m1 = T::from_w(coeffs.m1);
        let m2 = T::from_w(coeffs.m2);
        let fs = T::from_w(coeffs.fs);
        WideIIR2Coefficients {
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

    pub fn empty() -> WideIIR2Coefficients<T> {
        WideIIR2Coefficients {
            a: T::ZERO,
            g: T::ZERO,
            gpow2: T::ZERO,
            k: T::ZERO,
            a1: T::ZERO,
            a2: T::ZERO,
            a3: T::ZERO,
            m0: T::ZERO,
            m1: T::ZERO,
            m2: T::ZERO,
            fs: T::ZERO,
        }
    }

    pub fn empty_cascade() -> [WideIIR2Coefficients<T>; MAX_POLE_COUNT] {
        [WideIIR2Coefficients::<T>::empty(); MAX_POLE_COUNT]
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct WideIIR2<T: WIDE> {
    ic1eq: T,
    ic2eq: T,
    pub coeffs: WideIIR2Coefficients<T>,
}

impl<T: WIDE> WideIIR2<T> {
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: WideIIR2Coefficients<T>) -> Self {
        WideIIR2 {
            ic1eq: WIDE::ZERO,
            ic2eq: WIDE::ZERO,
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: T) -> T {
        let v3 = input - self.ic2eq;
        let v1 = self.coeffs.a1 * self.ic1eq + self.coeffs.a2 * v3;
        let v2 = self.ic2eq + self.coeffs.a2 * self.ic1eq + self.coeffs.a3 * v3;
        self.ic1eq = T::N2 * v1 - self.ic1eq;
        self.ic2eq = T::N2 * v2 - self.ic2eq;

        self.coeffs.m0 * input + self.coeffs.m1 * v1 + self.coeffs.m2 * v2
    }

    pub fn process_partial(&mut self, input: T) -> (T, T) {
        let v3 = input - self.ic2eq;
        let v1 = self.coeffs.a1 * self.ic1eq + self.coeffs.a2 * v3;
        let v2 = self.ic2eq + self.coeffs.a2 * self.ic1eq + self.coeffs.a3 * v3;
        self.ic1eq = T::N2 * v1 - self.ic1eq;
        self.ic2eq = T::N2 * v2 - self.ic2eq;

        (v1, v2)
    }

    pub fn update_coefficients(&mut self, new_coefficients: WideIIR2Coefficients<T>) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
mod tests {

    use wide::f64x4;

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
        let coeffs = WideIIR2Coefficients::from(coeffs);

        let mut filter_left = WideIIR2::new(coeffs);

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
