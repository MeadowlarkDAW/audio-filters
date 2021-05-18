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

    pub fn lowpass(f0: T, fs: T) -> IIR1Coefficients<T> {
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

    pub fn highpass(f0: T, fs: T) -> IIR1Coefficients<T> {
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

    pub fn allpass(f0: T, fs: T) -> IIR1Coefficients<T> {
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

use wide::f64x4;

#[derive(Copy, Clone, Debug)]
pub struct WideIIR1Coefficients {
    pub a: f64x4,
    pub g: f64x4,
    pub a1: f64x4,
    pub m0: f64x4,
    pub m1: f64x4,
    pub fs: f64x4,
}

impl WideIIR1Coefficients {
    pub fn from<T: FP>(coeffs: IIR1Coefficients<T>) -> WideIIR1Coefficients {
        let a = f64x4::splat(Into::<f64>::into(coeffs.a));
        let g = f64x4::splat(Into::<f64>::into(coeffs.g));
        let a1 = f64x4::splat(Into::<f64>::into(coeffs.a1));
        let m0 = f64x4::splat(Into::<f64>::into(coeffs.m0));
        let m1 = f64x4::splat(Into::<f64>::into(coeffs.m1));
        let fs = f64x4::splat(Into::<f64>::into(coeffs.fs));
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
pub struct WideIIR1 {
    ic1eq: f64x4,
    pub coeffs: WideIIR1Coefficients,
}

impl WideIIR1 {
    pub fn new(coefficients: WideIIR1Coefficients) -> Self {
        WideIIR1 {
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

    pub fn update_coefficients(&mut self, new_coefficients: WideIIR1Coefficients) {
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

        let coeffs = IIR1Coefficients::lowpass(f0, fs);
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
