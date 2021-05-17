use core::ops::{Add, Mul, Sub};
use num_complex::Complex;
use num_traits::{Float, FloatConst, One, Zero};

use crate::units::ZSample;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Errors {
    OutsideNyquist,
    NegativeQ,
    NegativeFrequency,
}

#[derive(Clone, Copy, Debug)]
pub enum IIR2Type {
    LowPass,
    HighPass,
    BandPass,
    Notch,
    AllPass,
    LowShelf,
    HighShelf,
    Bell,
}

#[derive(Copy, Clone, Debug)]
pub struct IIR2Coefficients<T> {
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

impl<T> IIR2Coefficients<T>
where
    T: Float,
    T: Zero,
    T: One,
    T: FloatConst,
    T: Into<Complex<T>>,
    f32: Into<T>,
    T: Add<Complex<T>, Output = Complex<T>>,
    T: Mul<Complex<T>, Output = Complex<T>>,
    T: Sub<Complex<T>, Output = Complex<T>>,
{
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

    /// Creates a SVF from a set of filter coefficients
    pub fn from_params(
        filter: IIR2Type,
        fs: T,
        f0: T,
        q_value: T,
        db_gain: T,
    ) -> Result<IIR2Coefficients<T>, Errors> {
        let two: T = 2.0.into();
        let ten: T = 10.0.into();
        let forty: T = 40.0.into();

        if two * f0 > fs {
            return Err(Errors::OutsideNyquist);
        }

        if q_value < T::zero() {
            return Err(Errors::NegativeQ);
        }

        let a;
        let g;
        let k;
        let a1;
        let a2;
        let a3;
        let m0;
        let m1;
        let m2;

        match filter {
            IIR2Type::LowPass
            | IIR2Type::HighPass
            | IIR2Type::AllPass
            | IIR2Type::BandPass
            | IIR2Type::Notch => {
                a = T::one();
                g = (T::PI() * f0 / fs).tan();
                k = T::one() / q_value;
                a1 = T::one() / (T::one() + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
            IIR2Type::LowShelf => {
                a = (ten).powf(db_gain / forty);
                g = (T::PI() * f0 / fs).tan() / a.sqrt();
                k = T::one() / q_value;
                a1 = T::one() / (T::one() + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
            IIR2Type::HighShelf => {
                a = ten.powf(db_gain / forty);
                g = (T::PI() * f0 / fs).tan() * a.sqrt();
                k = T::one() / q_value;
                a1 = T::one() / (T::one() + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
            IIR2Type::Bell => {
                a = ten.powf(db_gain / forty);
                g = (T::PI() * f0 / fs).tan();
                k = T::one() / (q_value * a);
                a1 = T::one() / (T::one() + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
        };

        match filter {
            IIR2Type::LowPass => {
                m0 = T::zero();
                m1 = T::zero();
                m2 = T::one();
            }
            IIR2Type::HighPass => {
                m0 = T::one();
                m1 = -k;
                m2 = -T::one();
            }
            IIR2Type::BandPass => {
                m0 = T::zero();
                m1 = T::one();
                m2 = T::zero();
            }
            IIR2Type::Notch => {
                m0 = T::one();
                m1 = -k;
                m2 = T::zero();
            }
            IIR2Type::AllPass => {
                m0 = T::one();
                m1 = -two * k;
                m2 = T::zero();
            }
            IIR2Type::LowShelf => {
                m0 = T::one();
                m1 = k * (a - T::one());
                m2 = a * a - T::one();
            }
            IIR2Type::HighShelf => {
                m0 = a * a;
                m1 = k * (T::one() - a) * a;
                m2 = T::one() - a * a;
            }
            IIR2Type::Bell => {
                m0 = T::one();
                m1 = k * (a * a - T::one());
                m2 = T::zero();
            }
        };
        Ok(IIR2Coefficients {
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
        })
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct IIR2<T> {
    ic1eq: T,
    ic2eq: T,
    pub coeffs: IIR2Coefficients<T>,
}

impl<T> IIR2<T>
where
    T: Float,
    T: Zero,
    f32: Into<T>,
{
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: IIR2Coefficients<T>) -> Self {
        IIR2 {
            ic1eq: T::zero(),
            ic2eq: T::zero(),
            coeffs: coefficients,
        }
    }

    pub fn run(&mut self, input: T) -> T {
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
