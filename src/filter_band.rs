use core::fmt;

use crate::units::FP;
use num_complex::Complex;

use num_traits::NumCast;

use crate::{
    first_order_iir::{IIR1Coefficients, IIR1},
    second_order_iir::{IIR2Coefficients, IIR2},
    units::{butterworth_cascade_q, Units, ZSample},
    MAX_POLE_COUNT,
};

#[derive(PartialEq, Debug, Clone, Copy)]
pub enum BandType {
    Bell,
    LowPass,
    HighPass,
    LowShelf,
    HighShelf,
    Notch,
    BandPass,
    AllPass,
}

impl BandType {
    pub fn from_u8(value: u8) -> BandType {
        match value {
            1 => BandType::Bell,
            2 => BandType::LowPass,
            3 => BandType::HighPass,
            4 => BandType::LowShelf,
            5 => BandType::HighShelf,
            6 => BandType::Notch,
            7 => BandType::BandPass,
            8 => BandType::AllPass,
            _ => BandType::LowPass,
        }
    }
}

impl fmt::Display for BandType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct IIRCoefficientsSet<T: FP> {
    pub iir1: IIR1Coefficients<T>,
    pub iir2: [IIR2Coefficients<T>; MAX_POLE_COUNT],
    pub iir1_enabled: bool,
}

impl<T: FP> IIRCoefficientsSet<T> {
    pub fn new(sample_rate: T) -> IIRCoefficientsSet<T> {
        let iir2_coeffs = IIR2Coefficients::bell(sample_rate, T::one(), 1000.0.into(), T::zero());
        let iir1_coeffs = IIR1Coefficients::<T>::lowpass(1000.0.into(), sample_rate);
        IIRCoefficientsSet {
            iir1_enabled: false,
            iir1: iir1_coeffs,
            iir2: [iir2_coeffs; MAX_POLE_COUNT],
        }
    }

    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        let mut y = if self.iir1_enabled {
            self.iir1.get_bode_sample(z.z) * self.iir2[0].get_bode_sample(z)
        } else {
            self.iir2[0].get_bode_sample(z)
        };
        for band_a in self.iir2.iter().skip(1) {
            y = y * band_a.get_bode_sample(z);
        }
        y
    }
}

#[derive(Copy, Clone, Debug)]
pub struct FilterBand<T: FP> {
    iir1: IIR1<T>,
    iir2: [IIR2<T>; MAX_POLE_COUNT],
    pub coeffs: IIRCoefficientsSet<T>,
    kind: BandType,
    freq: T,
    gain: T,
    bw_value: T,
    slope: T,
    u_slope: u8,
    iir1_enabled: bool,
    start_pole: usize,
    iir2_cascade_count: usize,
    sample_rate: T,
}

impl<T: FP> FilterBand<T> {
    pub fn new(sample_rate: T) -> FilterBand<T> {
        let coeffs = IIRCoefficientsSet::<T>::new(sample_rate);
        FilterBand {
            iir1: IIR1::<T>::new(coeffs.iir1),
            iir2: [IIR2::<T>::new(coeffs.iir2[0]); MAX_POLE_COUNT],
            coeffs,
            kind: BandType::Bell,
            freq: T::zero(),
            gain: T::zero(),
            bw_value: T::zero(),
            slope: 2.0.into(),
            u_slope: 2,
            iir1_enabled: false,
            start_pole: 0,
            iir2_cascade_count: 1,
            sample_rate: 48000.0.into(),
        }
    }

    pub fn get_bode_sample(&self, z: ZSample<T>) -> Complex<T> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        self.coeffs.get_bode_sample(z)
    }

    pub fn process(&mut self, x: T) -> T {
        let mut x = x;

        if self.iir1_enabled {
            x = self.iir1.run(x);
        }
        if self.u_slope > 1 {
            if self.kind == BandType::Bell || self.kind == BandType::Notch {
                x = self.iir2[0].run(x);
            } else {
                for i in self.start_pole..self.iir2_cascade_count {
                    x = self.iir2[i].run(x);
                }
            }
        }
        x
    }

    pub fn update_coeffs(&mut self, coeffs: IIRCoefficientsSet<T>) {
        for (filter, coeff) in self.iir2.iter_mut().zip(coeffs.iir2.iter()) {
            filter.update_coefficients(*coeff)
        }
        self.iir1.update_coefficients(coeffs.iir1);
    }

    pub fn mimic_band(&mut self, band: &FilterBand<T>) {
        self.coeffs = band.coeffs;
        self.kind = band.kind;
        self.freq = band.freq;
        self.gain = band.gain;
        self.bw_value = band.bw_value;
        self.slope = band.slope;
        self.u_slope = band.u_slope;
        self.iir1_enabled = band.iir1_enabled;
        self.start_pole = band.start_pole;
        self.iir2_cascade_count = band.iir2_cascade_count;
        self.sample_rate = band.sample_rate;
        self.update_coeffs(self.coeffs);
    }

    pub fn update(
        &mut self,
        kind: BandType,
        in_freq: T,
        in_gain: T,
        in_bw_value: T,
        slope: T,
        fs: T,
    ) {
        if kind == self.kind
            && in_freq == self.freq
            && in_gain == self.gain
            && in_bw_value == self.bw_value
            && slope == self.slope
            && fs == self.sample_rate
        {
            return;
        }

        self.kind = kind;
        self.freq = in_freq;
        self.gain = in_gain;
        self.bw_value = in_bw_value;
        self.slope = slope;
        self.sample_rate = fs;

        let f0 = self.freq;
        let gain = self.gain;
        let bw_value = self.bw_value;

        self.u_slope = NumCast::from(slope).unwrap();

        let odd_order = self.u_slope & 1 == 1;

        self.iir1_enabled = odd_order
            && match kind {
                BandType::LowPass => true,
                BandType::HighPass => true,
                BandType::LowShelf => true,
                BandType::HighShelf => true,
                BandType::AllPass => true,
                BandType::BandPass => true,
                _ => false,
            };

        self.coeffs.iir1_enabled = self.iir1_enabled;

        self.start_pole = if self.iir1_enabled { 1 } else { 0 } as usize;
        self.iir2_cascade_count = ((self.u_slope as usize + self.start_pole) / 2) as usize;

        let q_value = bw_value.bw_to_q(f0, fs);
        let q_offset = q_value * T::FRAC_1_SQRT_2(); //butterworth Q

        let mut partial_gain = gain / self.u_slope.into();

        if self.iir1_enabled {
            self.coeffs.iir1 = match kind {
                BandType::LowPass => IIR1Coefficients::lowpass(f0, fs),
                BandType::HighPass => IIR1Coefficients::highpass(f0, fs),
                BandType::LowShelf => IIR1Coefficients::lowshelf(f0, partial_gain, fs),
                BandType::HighShelf => IIR1Coefficients::highshelf(f0, partial_gain, fs),
                BandType::AllPass => IIR1Coefficients::allpass(f0, fs),
                _ => IIR1Coefficients::lowpass(f0, fs),
            };
        }

        if self.u_slope <= 1 {
            self.update_coeffs(self.coeffs);
            return;
        }

        let two: T = 2.0.into();
        partial_gain = partial_gain * two;

        match self.kind {
            BandType::Bell => {
                self.coeffs.iir2[0] = IIR2Coefficients::bell(f0, q_value, gain, fs);
            }
            BandType::LowPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
                    self.coeffs.iir2[i] = IIR2Coefficients::lowpass(f0, q_value * q_offset, fs);
                }
            }
            BandType::HighPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
                    self.coeffs.iir2[i] = IIR2Coefficients::highpass(f0, q_value * q_offset, fs);
                }
            }
            BandType::LowShelf => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
                    self.coeffs.iir2[i] =
                        IIR2Coefficients::lowshelf(f0, q_value * q_offset, partial_gain, fs);
                }
            }
            BandType::HighShelf => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
                    self.coeffs.iir2[i] =
                        IIR2Coefficients::highshelf(f0, q_value * q_offset, partial_gain, fs);
                }
            }
            BandType::Notch => {
                self.coeffs.iir2[0] = IIR2Coefficients::notch(f0, q_value * q_offset, fs);
            }
            BandType::BandPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
                    self.coeffs.iir2[i] = IIR2Coefficients::bandpass(f0, q_value * q_offset, fs);
                }
            }
            BandType::AllPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value: T = butterworth_cascade_q(self.u_slope, i as u8);
                    self.coeffs.iir2[i] = IIR2Coefficients::allpass(f0, q_value * q_offset, fs);
                }
            }
        }
        self.update_coeffs(self.coeffs);
    }
}
