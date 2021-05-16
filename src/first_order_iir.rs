use core::f64::consts::{PI, TAU};

use libm::*;
use num_complex::Complex;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Errors {
    OutsideNyquist,
    NegativeQ,
    NegativeFrequency,
}

pub fn get_z(f_hz: f64, fs: f64) -> Complex<f64> {
    let z = -TAU * f_hz / fs;
    cos(z) + sin(z) * Complex::new(0.0, 1.0)
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
pub struct IIR1Coefficients<T> {
    pub a: T,
    pub g: T,
    pub a1: T,
    pub m0: T,
    pub m1: T,
    pub fs: T,
}

impl IIR1Coefficients<f64> {
    pub fn get_bode_sample(self, z: Complex<f64>) -> Complex<f64> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let denominator = self.g + z * (self.g - 1.0) + 1.0;

        let y = self.m0 + (self.m1 * self.g * (z + 1.0)) / denominator;

        y
    }

    /// Creates a SVF from a set of filter coefficients
    pub fn from_params(
        filter: IIR1Type<f64>,
        fs: f64,
        f0: f64,
    ) -> Result<IIR1Coefficients<f64>, Errors> {
        if 2.0 * f0 > fs {
            return Err(Errors::OutsideNyquist);
        }

        let a;
        let g;
        let a1;
        let m0;
        let m1;

        match filter {
            IIR1Type::LowPass | IIR1Type::HighPass | IIR1Type::AllPass => {
                a = 1.0;
                g = (PI * f0 / fs).tan();
                a1 = g / (1.0 + g);
            }
            IIR1Type::LowShelf(db_gain) => {
                a = 10.0f64.powf(db_gain / 20.0);
                g = (PI * f0 / fs).tan() / (a).sqrt();
                a1 = g / (1.0 + g);
            }
            IIR1Type::HighShelf(db_gain) => {
                a = 10.0f64.powf(db_gain / 20.0);
                g = (PI * f0 / fs).tan() * (a).sqrt();
                a1 = g / (1.0 + g);
            }
        };

        match filter {
            IIR1Type::LowPass => {
                m0 = 0.0;
                m1 = 1.0;
            }
            IIR1Type::HighPass => {
                m0 = 1.0;
                m1 = -1.0;
            }
            IIR1Type::AllPass => {
                m0 = 1.0;
                m1 = -2.0;
            }
            IIR1Type::LowShelf(_) => {
                m0 = 1.0;
                m1 = a - 1.0;
            }
            IIR1Type::HighShelf(_) => {
                m0 = a;
                m1 = 1.0 - a;
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
pub struct IIR1<T> {
    ic1eq: T,
    pub coeffs: IIR1Coefficients<T>,
}

impl IIR1<f64> {
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: IIR1Coefficients<f64>) -> Self {
        IIR1 {
            ic1eq: 0.0,
            coeffs: coefficients,
        }
    }

    pub fn run(&mut self, input: f64) -> f64 {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input + self.coeffs.m1 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR1Coefficients<f64>) {
        self.coeffs = new_coefficients;
    }
}
