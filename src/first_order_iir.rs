use num_complex::Complex;

use crate::units::FP;

pub fn get_z<T: FP>(f_hz: T, fs: T) -> Complex<T> {
    let z = -T::TAU() * f_hz / fs;
    z.cos() + z.sin() * Complex::<T>::new(T::zero(), T::one())
}

#[derive(Copy, Clone, Debug)]
pub struct IIR1Coefficients<T: FP> {
    pub a: T,
    pub g: T,
    pub a1: T,
    pub m0: T,
    pub m1: T,
    pub fs: T,
}

impl<T: FP> IIR1Coefficients<T> {
    pub fn get_bode_sample(self, z: Complex<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let denominator = self.g + z * (self.g - T::one()) + T::one();

        let y = self.m0 + (self.m1 * self.g * (z + T::one())) / denominator;

        y
    }

    pub fn empty() -> IIR1Coefficients<T> {
        IIR1Coefficients {
            a: T::zero(),
            g: T::zero(),
            a1: T::zero(),
            m0: T::zero(),
            m1: T::zero(),
            fs: T::zero(),
        }
    }

    pub fn lowpass(f0: T, _db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let a1 = g / (T::one() + g);
        let m0 = T::zero();
        let m1 = T::one();
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn highpass(f0: T, _db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let a1 = g / (T::one() + g);
        let m0 = T::one();
        let m1 = -T::one();
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn allpass(f0: T, _db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = T::one();
        let g = (T::PI() * f0 / fs).tan();
        let a1 = g / (T::one() + g);
        let m0 = T::one();
        let m1 = -Into::<T>::into(2.0);
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn lowshelf(f0: T, db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = Into::<T>::into(10.0).powf(db_gain / Into::<T>::into(20.0));
        let g = (T::PI() * f0 / fs).tan() / (a).sqrt();
        let a1 = g / (T::one() + g);
        let m0 = T::one();
        let m1 = a - T::one();
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn highshelf(f0: T, db_gain: T, fs: T) -> IIR1Coefficients<T> {
        let f0 = f0.min(fs * Into::<T>::into(0.5));
        let a = Into::<T>::into(10.0).powf(db_gain / Into::<T>::into(20.0));
        let g = (T::PI() * f0 / fs).tan() * (a).sqrt();
        let a1 = g / (T::one() + g);
        let m0 = a;
        let m1 = T::one() - a;
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
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
            ic1eq: T::zero(),
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: T) -> T {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input + self.coeffs.m1 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR1Coefficients<T>) {
        self.coeffs = new_coefficients;
    }
}
