use num_complex::Complex;

use crate::{units::ZSample, MAX_POLE_COUNT};

use crate::units::FP;

#[derive(Copy, Clone, Debug)]
pub struct IIR2Coefficients<T: FP> {
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

impl<T: FP> IIR2Coefficients<T> {
    pub fn get_bode_sample(self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let two: T = 2.0.into();

        let denominator = (self.gpow2 + self.g * self.k + T::one())
            + two * (self.gpow2 - T::one()) * z.pow1
            + (self.gpow2 - self.g * self.k + T::one()) * z.pow2;

        let y = self.m0
            + (self.m1 * self.g * (T::one() - z.pow2)
                + self.m2 * self.gpow2 * (T::one() + two * z.pow1 + z.pow2))
                / denominator;

        y
    }

    pub fn empty() -> IIR2Coefficients<T> {
        IIR2Coefficients {
            a: T::zero(),
            g: T::zero(),
            gpow2: T::zero(),
            k: T::zero(),
            a1: T::zero(),
            a2: T::zero(),
            a3: T::zero(),
            m0: T::zero(),
            m1: T::zero(),
            m2: T::zero(),
            fs: T::zero(),
        }
    }

    pub fn empty_cascade() -> [IIR2Coefficients<T>; MAX_POLE_COUNT] {
        [IIR2Coefficients::empty(); MAX_POLE_COUNT]
    }

    pub fn lowpass(f0: T, q_value: T, _db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::zero();
        let m1 = T::zero();
        let m2 = T::one();
        IIR2Coefficients {
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
    pub fn highpass(f0: T, q_value: T, _db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::one();
        let m1 = -k;
        let m2 = -T::one();
        IIR2Coefficients {
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
    pub fn bandpass(f0: T, q_value: T, _db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::zero();
        let m1 = T::one();
        let m2 = T::zero();
        IIR2Coefficients {
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
    pub fn notch(f0: T, q_value: T, _db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::one();
        let m1 = -k;
        let m2 = T::zero();
        IIR2Coefficients {
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
    pub fn allpass(f0: T, q_value: T, _db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::one();
        let m1 = -Into::<T>::into(2.0) * k;
        let m2 = T::zero();
        IIR2Coefficients {
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
    pub fn lowshelf(f0: T, q_value: T, db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = Into::<T>::into(10.0).powf(db_gain / Into::<T>::into(40.0));
        let g = (T::PI() * f0 / fs).tan() / a.sqrt();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::one();
        let m1 = k * (a - T::one());
        let m2 = a * a - T::one();
        IIR2Coefficients {
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
    pub fn highshelf(f0: T, q_value: T, db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = Into::<T>::into(10.0).powf(db_gain / Into::<T>::into(40.0));
        let g = (T::PI() * f0 / fs).tan() * a.sqrt();
        let k = T::one() / q_value;
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = a * a;
        let m1 = k * (T::one() - a) * a;
        let m2 = T::one() - a * a;
        IIR2Coefficients {
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
    pub fn bell(f0: T, q_value: T, db_gain: T, fs: T) -> IIR2Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = Into::<T>::into(10.0).powf(db_gain / Into::<T>::into(40.0));
        let g = (T::PI() * f0 / fs).tan();
        let k = T::one() / (q_value * a);
        let a1 = T::one() / (T::one() + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = T::one();
        let m1 = k * (a * a - T::one());
        let m2 = T::zero();
        IIR2Coefficients {
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
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct IIR2<T: FP> {
    ic1eq: T,
    ic2eq: T,
    pub coeffs: IIR2Coefficients<T>,
}

impl<T: FP> IIR2<T> {
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: IIR2Coefficients<T>) -> Self {
        IIR2 {
            ic1eq: T::zero(),
            ic2eq: T::zero(),
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: T) -> T {
        let two: T = 2.0.into();
        let v3 = input - self.ic2eq;
        let v1 = self.coeffs.a1 * self.ic1eq + self.coeffs.a2 * v3;
        let v2 = self.ic2eq + self.coeffs.a2 * self.ic1eq + self.coeffs.a3 * v3;
        self.ic1eq = two * v1 - self.ic1eq;
        self.ic2eq = two * v2 - self.ic2eq;

        self.coeffs.m0 * input + self.coeffs.m1 * v1 + self.coeffs.m2 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR2Coefficients<T>) {
        self.coeffs = new_coefficients;
    }
}
