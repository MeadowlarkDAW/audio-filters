use num_complex::Complex;

use crate::units::FP;

pub fn get_z<T: FP>(f_hz: T, fs: T) -> Complex<T> {
    let z = -T::TAU() * f_hz / fs;
    z.cos() + z.sin() * Complex::<T>::new(T::N0, T::N1)
}

#[derive(Copy, Clone, Debug)]
pub struct IIR1Coefficients<T: FP> {
    pub a: T,
    pub g: T,
    pub a1: T,
    pub m0: T,
    pub m1: T,
}

impl<T: FP> IIR1Coefficients<T> {
    pub fn get_bode_sample(self, z: Complex<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let denominator = self.g + z * (self.g - T::N1) + T::N1;

        let y = self.m0 + (self.m1 * self.g * (z + T::N1)) / denominator;

        y
    }

    //TODO make const once possible
    pub fn empty() -> IIR1Coefficients<T> {
        IIR1Coefficients {
            a: T::N0,
            g: T::N0,
            a1: T::N0,
            m0: T::N0,
            m1: T::N0,
        }
    }

    pub fn lowpass(f0: T, _db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * T::N0_5);
        let a = T::N1;
        let g = (T::PI() * f0 / fs).tan();
        let a1 = g / (T::N1 + g);
        let m0 = T::N0;
        let m1 = T::N1;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    pub fn highpass(f0: T, _db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * T::N0_5);
        let a = T::N1;
        let g = (T::PI() * f0 / fs).tan();
        let a1 = g / (T::N1 + g);
        let m0 = T::N1;
        let m1 = -T::N1;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    pub fn allpass(f0: T, _db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * T::N0_5);
        let a = T::N1;
        let g = (T::PI() * f0 / fs).tan();
        let a1 = g / (T::N1 + g);
        let m0 = T::N1;
        let m1 = -T::N2;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    pub fn lowshelf(f0: T, db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * T::N0_5);
        let a = T::N10.powf(db_gain / T::N20);
        let g = (T::PI() * f0 / fs).tan() / (a).sqrt();
        let a1 = g / (T::N1 + g);
        let m0 = T::N1;
        let m1 = a - T::N1;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    pub fn highshelf(f0: T, db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * T::N0_5);
        let a = T::N10.powf(db_gain / T::N20);
        let g = (T::PI() * f0 / fs).tan() * (a).sqrt();
        let a1 = g / (T::N1 + g);
        let m0 = a;
        let m1 = T::N1 - a;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct IIR1<T: FP> {
    ic1eq: T,
    pub coeffs: IIR1Coefficients<T>,
}

impl<T: FP> IIR1<T> {
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: IIR1Coefficients<T>) -> Self {
        IIR1 {
            ic1eq: T::N0,
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

    pub fn update_coefficients(&mut self, new_coefficients: IIR1Coefficients<T>) {
        self.coeffs = new_coefficients;
    }
}
