use core::{f64::consts::FRAC_1_SQRT_2, fmt};

use num_complex::Complex;

use crate::{
    first_order_iir::{IIR1Coefficients, IIR1Type, IIR1},
    second_order_iir::{IIR2Coefficients, IIR2Type, IIR2},
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
    pub fn from_u32(value: u32) -> BandType {
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
pub struct IIRCoefficientsSet {
    pub iir1: IIR1Coefficients<f64>,
    pub iir2: [IIR2Coefficients<f64>; MAX_POLE_COUNT],
    pub iir1_enabled: bool,
}

impl IIRCoefficientsSet {
    pub fn new(sample_rate: f64) -> IIRCoefficientsSet {
        let iir2_coeffs =
            IIR2Coefficients::<f64>::from_params(IIR2Type::Bell(0.0f64), sample_rate, 1000.0, 1.0)
                .unwrap();
        let iir1_coeffs =
            IIR1Coefficients::<f64>::from_params(IIR1Type::LowPass, sample_rate, 1000.0).unwrap();
        IIRCoefficientsSet {
            iir1_enabled: false,
            iir1: iir1_coeffs,
            iir2: [iir2_coeffs; MAX_POLE_COUNT],
        }
    }

    pub fn get_bode_sample(&self, z: ZSample<f64>) -> Complex<f64> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        let mut y = if self.iir1_enabled {
            self.iir1.get_bode_sample(z.z) * self.iir2[0].get_bode_sample(z)
        } else {
            self.iir2[0].get_bode_sample(z)
        };
        for band_a in self.iir2.iter().skip(1) {
            y *= band_a.get_bode_sample(z);
        }
        y
    }
}

#[derive(Copy, Clone, Debug)]
pub struct FilterBand {
    iir1: IIR1<f64>,
    iir2: [IIR2<f64>; MAX_POLE_COUNT],
    pub coeffs: IIRCoefficientsSet,
    kind: BandType,
    freq: f64,
    gain: f64,
    bw_value: f64,
    slope: f64,
    u_slope: u32,
    iir1_enabled: bool,
    start_pole: usize,
    iir2_cascade_count: usize,
    sample_rate: f64,
}

impl FilterBand {
    pub fn new(sample_rate: f64) -> FilterBand {
        let coeffs = IIRCoefficientsSet::new(sample_rate);
        FilterBand {
            iir1: IIR1::<f64>::new(coeffs.iir1),
            iir2: [IIR2::<f64>::new(coeffs.iir2[0]); MAX_POLE_COUNT],
            coeffs,
            kind: BandType::Bell,
            freq: 0.0,
            gain: 0.0,
            bw_value: 0.0,
            slope: 2.0,
            u_slope: 2,
            iir1_enabled: false,
            start_pole: 0,
            iir2_cascade_count: 1,
            sample_rate: 48000.0,
        }
    }

    pub fn get_bode_sample(&self, z: ZSample<f64>) -> Complex<f64> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.
        self.coeffs.get_bode_sample(z)
    }

    pub fn process(&mut self, x: f64) -> f64 {
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

    pub fn update_coeffs(&mut self, coeffs: IIRCoefficientsSet) {
        for (filter, coeff) in self.iir2.iter_mut().zip(coeffs.iir2.iter()) {
            filter.update_coefficients(*coeff)
        }
        self.iir1.update_coefficients(coeffs.iir1);
    }

    pub fn mimic_band(&mut self, band: &FilterBand) {
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
        in_freq: f64,
        in_gain: f64,
        in_bw_value: f64,
        slope: f64,
        sample_rate: f64,
    ) {
        if kind == self.kind
            && in_freq == self.freq
            && in_gain == self.gain
            && in_bw_value == self.bw_value
            && slope == self.slope
            && sample_rate == self.sample_rate
        {
            return;
        }

        self.kind = kind;
        self.freq = in_freq;
        self.gain = in_gain;
        self.bw_value = in_bw_value;
        self.slope = slope;
        self.sample_rate = sample_rate;

        let freq = self.freq;
        let gain = self.gain;
        let bw_value = self.bw_value;

        self.u_slope = slope as u32;

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

        let q_value = bw_value.bw_to_q(freq, sample_rate);
        let q_offset = q_value * FRAC_1_SQRT_2; //butterworth Q

        let mut partial_gain = gain / self.u_slope as f64;

        if self.iir1_enabled {
            self.coeffs.iir1 = IIR1Coefficients::<f64>::from_params(
                match kind {
                    BandType::LowPass => IIR1Type::LowPass,
                    BandType::HighPass => IIR1Type::HighPass,
                    BandType::LowShelf => IIR1Type::LowShelf(partial_gain),
                    BandType::HighShelf => IIR1Type::HighShelf(partial_gain),
                    BandType::AllPass => IIR1Type::AllPass,
                    BandType::BandPass => IIR1Type::HighPass,
                    _ => IIR1Type::LowPass,
                },
                sample_rate,
                freq,
            )
            .unwrap();
        }

        if self.u_slope <= 1 {
            self.update_coeffs(self.coeffs);
            return;
        }

        partial_gain *= 2.0;

        match self.kind {
            BandType::Bell => {
                self.coeffs.iir2[0] = IIR2Coefficients::<f64>::from_params(
                    IIR2Type::Bell(gain),
                    sample_rate,
                    freq,
                    q_value,
                )
                .unwrap();
            }
            BandType::LowPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value = butterworth_cascade_q(self.u_slope, i as u32);
                    self.coeffs.iir2[i] = IIR2Coefficients::<f64>::from_params(
                        IIR2Type::LowPass,
                        sample_rate,
                        freq,
                        q_value * q_offset,
                    )
                    .unwrap();
                }
            }
            BandType::HighPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value = butterworth_cascade_q(self.u_slope, i as u32);
                    self.coeffs.iir2[i] = IIR2Coefficients::<f64>::from_params(
                        IIR2Type::HighPass,
                        sample_rate,
                        freq,
                        q_value * q_offset,
                    )
                    .unwrap();
                }
            }
            BandType::LowShelf => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value = butterworth_cascade_q(self.u_slope, i as u32);
                    self.coeffs.iir2[i] = IIR2Coefficients::<f64>::from_params(
                        IIR2Type::LowShelf(partial_gain),
                        sample_rate,
                        freq,
                        q_value * q_offset,
                    )
                    .unwrap();
                }
            }
            BandType::HighShelf => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value = butterworth_cascade_q(self.u_slope, i as u32);
                    self.coeffs.iir2[i] = IIR2Coefficients::<f64>::from_params(
                        IIR2Type::HighShelf(partial_gain),
                        sample_rate,
                        freq,
                        q_value * q_offset,
                    )
                    .unwrap();
                }
            }
            BandType::Notch => {
                self.coeffs.iir2[0] = IIR2Coefficients::<f64>::from_params(
                    IIR2Type::Notch,
                    sample_rate,
                    freq,
                    q_value,
                )
                .unwrap();
            }
            BandType::BandPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value = butterworth_cascade_q(self.u_slope, i as u32);
                    self.coeffs.iir2[i] = IIR2Coefficients::<f64>::from_params(
                        IIR2Type::BandPass,
                        sample_rate,
                        freq,
                        q_value * q_offset,
                    )
                    .unwrap();
                }
            }
            BandType::AllPass => {
                for i in self.start_pole..self.iir2_cascade_count {
                    let q_value = butterworth_cascade_q(self.u_slope, i as u32);
                    self.coeffs.iir2[i] = IIR2Coefficients::<f64>::from_params(
                        IIR2Type::AllPass,
                        sample_rate,
                        freq,
                        q_value * q_offset,
                    )
                    .unwrap();
                }
            }
        }
        self.update_coeffs(self.coeffs);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_band() {
        let mut band = FilterBand::new(48000.0);
        band.update(BandType::LowPass, 1000.0, 0.0, 1.0, 4.0, 48000.0);
        dbg!(band.iir1_enabled);
        dbg!(band.u_slope);
        dbg!(band.start_pole);
        dbg!(band.iir2_cascade_count);
    }
}
