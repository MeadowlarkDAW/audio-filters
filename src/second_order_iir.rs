use core::f64::consts::PI;

use libm::*;
use num_complex::Complex;

use crate::units::ZSample;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Errors {
    OutsideNyquist,
    NegativeQ,
    NegativeFrequency,
}

#[derive(Clone, Copy, Debug)]
pub enum IIR2Type<DBGain> {
    LowPass,
    HighPass,
    BandPass,
    Notch,
    AllPass,
    LowShelf(DBGain),
    HighShelf(DBGain),
    Bell(DBGain),
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

impl IIR2Coefficients<f64> {
    pub fn get_bode_sample(self, z: ZSample<f64>) -> Complex<f64> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let denominator = (self.gpow2 + self.g * self.k + 1.0)
            + 2.0 * (self.gpow2 - 1.0) * z.pow1
            + (self.gpow2 - self.g * self.k + 1.0) * z.pow2;

        let y = self.m0
            + (self.m1 * self.g * (1.0 - z.pow2)
                + self.m2 * self.gpow2 * (1.0 + 2.0 * z.pow1 + z.pow2))
                / denominator;

        y
    }

    /// Creates a SVF from a set of filter coefficients
    pub fn from_params(
        filter: IIR2Type<f64>,
        fs: f64,
        f0: f64,
        q_value: f64,
    ) -> Result<IIR2Coefficients<f64>, Errors> {
        if 2.0 * f0 > fs {
            return Err(Errors::OutsideNyquist);
        }

        if q_value < 0.0 {
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
                a = 1.0;
                g = tan(PI * f0 / fs);
                k = 1.0 / q_value;
                a1 = 1.0 / (1.0 + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
            IIR2Type::LowShelf(db_gain) => {
                a = pow(10.0f64, db_gain / 40.0);
                g = tan(PI * f0 / fs) / sqrt(a);
                k = 1.0 / q_value;
                a1 = 1.0 / (1.0 + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
            IIR2Type::HighShelf(db_gain) => {
                a = pow(10.0f64, db_gain / 40.0);
                g = tan(PI * f0 / fs) * sqrt(a);
                k = 1.0 / q_value;
                a1 = 1.0 / (1.0 + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
            IIR2Type::Bell(db_gain) => {
                a = pow(10.0f64, db_gain / 40.0);
                g = tan(PI * f0 / fs);
                k = 1.0 / (q_value * a);
                a1 = 1.0 / (1.0 + g * (g + k));
                a2 = g * a1;
                a3 = g * a2;
            }
        };

        match filter {
            IIR2Type::LowPass => {
                m0 = 0.0;
                m1 = 0.0;
                m2 = 1.0;
            }
            IIR2Type::HighPass => {
                m0 = 1.0;
                m1 = -k;
                m2 = -1.0;
            }
            IIR2Type::BandPass => {
                m0 = 0.0;
                m1 = 1.0;
                m2 = 0.0;
            }
            IIR2Type::Notch => {
                m0 = 1.0;
                m1 = -k;
                m2 = 0.0;
            }
            IIR2Type::AllPass => {
                m0 = 1.0;
                m1 = -2.0 * k;
                m2 = 0.0;
            }
            IIR2Type::LowShelf(_) => {
                m0 = 1.0;
                m1 = k * (a - 1.0);
                m2 = a * a - 1.0;
            }
            IIR2Type::HighShelf(_) => {
                m0 = a * a;
                m1 = k * (1.0 - a) * a;
                m2 = 1.0 - a * a;
            }
            IIR2Type::Bell(_) => {
                m0 = 1.0;
                m1 = k * (a * a - 1.0);
                m2 = 0.0;
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

impl IIR2<f64> {
    /// Creates a SVF from a set of filter coefficients
    pub fn new(coefficients: IIR2Coefficients<f64>) -> Self {
        IIR2 {
            ic1eq: 0.0,
            ic2eq: 0.0,
            coeffs: coefficients,
        }
    }

    pub fn run(&mut self, input: f64) -> f64 {
        let v3 = input - self.ic2eq;
        let v1 = self.coeffs.a1 * self.ic1eq + self.coeffs.a2 * v3;
        let v2 = self.ic2eq + self.coeffs.a2 * self.ic1eq + self.coeffs.a3 * v3;
        self.ic1eq = 2.0 * v1 - self.ic1eq;
        self.ic2eq = 2.0 * v2 - self.ic2eq;

        self.coeffs.m0 * input + self.coeffs.m1 * v1 + self.coeffs.m2 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR2Coefficients<f64>) {
        self.coeffs = new_coefficients;
    }
}
