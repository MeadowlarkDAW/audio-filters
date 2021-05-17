use num_complex::Complex;

use crate::units::FP;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Errors {
    OutsideNyquist,
    NegativeQ,
    NegativeFrequency,
}

pub fn get_z<T: FP>(f_hz: T, fs: T) -> Complex<T> {
    let z = -T::TAU() * f_hz / fs;
    z.cos() + z.sin() * Complex::<T>::new(T::zero(), T::one())
}

#[derive(Clone, Copy, Debug)]
pub enum IIR1Type<DBGain> {
    LowPass,
    HighPass,
    AllPass,
    LowShelf(DBGain),
    HighShelf(DBGain),
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

    /// Creates a SVF from a set of filter coefficients
    pub fn from_params(filter: IIR1Type<T>, fs: T, f0: T) -> Result<IIR1Coefficients<T>, Errors> {
        let two: T = 2.0.into();
        let ten: T = 10.0.into();
        let twenty: T = 20.0.into();

        if two * f0 > fs {
            return Err(Errors::OutsideNyquist);
        }

        let a;
        let g;
        let a1;
        let m0;
        let m1;

        match filter {
            IIR1Type::LowPass | IIR1Type::HighPass | IIR1Type::AllPass => {
                a = T::one();
                g = (T::PI() * f0 / fs).tan();
                a1 = g / (T::one() + g);
            }
            IIR1Type::LowShelf(db_gain) => {
                a = ten.powf(db_gain / twenty);
                g = (T::PI() * f0 / fs).tan() / (a).sqrt();
                a1 = g / (T::one() + g);
            }
            IIR1Type::HighShelf(db_gain) => {
                a = ten.powf(db_gain / twenty);
                g = (T::PI() * f0 / fs).tan() * (a).sqrt();
                a1 = g / (T::one() + g);
            }
        };

        match filter {
            IIR1Type::LowPass => {
                m0 = T::zero();
                m1 = T::one();
            }
            IIR1Type::HighPass => {
                m0 = T::one();
                m1 = -T::one();
            }
            IIR1Type::AllPass => {
                m0 = T::one();
                m1 = -two;
            }
            IIR1Type::LowShelf(_) => {
                m0 = T::one();
                m1 = a - T::one();
            }
            IIR1Type::HighShelf(_) => {
                m0 = a;
                m1 = T::one() - a;
            }
        };
        Ok(IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        })
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

    pub fn run(&mut self, input: T) -> T {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input + self.coeffs.m1 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR1Coefficients<T>) {
        self.coeffs = new_coefficients;
    }
}
