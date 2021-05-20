use crate::{
    filter_band::ProcessType,
    filter_band_wide::{
        WideF32FilterBand, WideF32FilterBandCoefficients, WideF64FilterBand,
        WideF64FilterBandCoefficients,
    },
};

use wide::f32x8;
use wide::f64x4;

#[derive(Copy, Clone, Debug)]
pub struct WideF64LinkwitzRileyCoefficients {
    pub coeffs: WideF64FilterBandCoefficients,
}

impl WideF64LinkwitzRileyCoefficients {
    //The resulting Linkwitz-Riley filter will have 2x to order of the input coefficients and 2x gain
    pub fn from(coeffs: WideF64FilterBandCoefficients) -> Self {
        WideF64LinkwitzRileyCoefficients { coeffs }
    }
}

#[derive(Copy, Clone)]
pub struct WideF64LinkwitzRileyBand {
    pub filter1: WideF64FilterBand,
    pub filter2: WideF64FilterBand,
    pub process: fn(&mut Self, f64x4) -> f64x4,
}

impl WideF64LinkwitzRileyBand {
    pub fn from(lw_coeffs: &WideF64LinkwitzRileyCoefficients) -> WideF64LinkwitzRileyBand {
        WideF64LinkwitzRileyBand {
            filter1: WideF64FilterBand::from(&lw_coeffs.coeffs),
            filter2: WideF64FilterBand::from(&lw_coeffs.coeffs),
            process: WideF64LinkwitzRileyBand::get_process(lw_coeffs.coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, x: f64x4) -> f64x4 {
        self.filter2
            .process_iir1_only(self.filter1.process_iir1_only(x))
    }

    pub fn process_iir2_only(&mut self, x: f64x4) -> f64x4 {
        self.filter2
            .process_iir2_only(self.filter1.process_iir2_only(x))
    }

    pub fn process_even_order_cascade(&mut self, x: f64x4) -> f64x4 {
        self.filter2
            .process_even_order_cascade(self.filter1.process_even_order_cascade(x))
    }

    pub fn process_odd_order_cascade(&mut self, x: f64x4) -> f64x4 {
        self.filter2
            .process_odd_order_cascade(self.filter1.process_odd_order_cascade(x))
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, f64x4) -> f64x4 {
        match process_type {
            ProcessType::ProcessIIR1Only => WideF64LinkwitzRileyBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => WideF64LinkwitzRileyBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => {
                WideF64LinkwitzRileyBand::process_even_order_cascade
            }
            ProcessType::ProcessOddOrderCascade => {
                WideF64LinkwitzRileyBand::process_odd_order_cascade
            }
        }
    }

    pub fn update(&mut self, lw_coeffs: &WideF64LinkwitzRileyCoefficients) {
        self.filter1.update(&lw_coeffs.coeffs);
        self.filter2.update(&lw_coeffs.coeffs);
        self.process = WideF64LinkwitzRileyBand::get_process(lw_coeffs.coeffs.process);
    }
}

//Copy of WideF64 (x4) just replacing with F32 (x8)
#[derive(Copy, Clone, Debug)]
pub struct WideF32LinkwitzRileyCoefficients {
    pub coeffs: WideF32FilterBandCoefficients,
}

impl WideF32LinkwitzRileyCoefficients {
    //The resulting Linkwitz-Riley filter will have 2x to order of the input coefficients and 2x gain
    pub fn from(coeffs: WideF32FilterBandCoefficients) -> Self {
        WideF32LinkwitzRileyCoefficients { coeffs }
    }
}

#[derive(Copy, Clone)]
pub struct WideF32LinkwitzRileyBand {
    pub filter1: WideF32FilterBand,
    pub filter2: WideF32FilterBand,
    pub process: fn(&mut Self, f32x8) -> f32x8,
}

impl WideF32LinkwitzRileyBand {
    pub fn from(lw_coeffs: &WideF32LinkwitzRileyCoefficients) -> WideF32LinkwitzRileyBand {
        WideF32LinkwitzRileyBand {
            filter1: WideF32FilterBand::from(&lw_coeffs.coeffs),
            filter2: WideF32FilterBand::from(&lw_coeffs.coeffs),
            process: WideF32LinkwitzRileyBand::get_process(lw_coeffs.coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, x: f32x8) -> f32x8 {
        self.filter2
            .process_iir1_only(self.filter1.process_iir1_only(x))
    }

    pub fn process_iir2_only(&mut self, x: f32x8) -> f32x8 {
        self.filter2
            .process_iir2_only(self.filter1.process_iir2_only(x))
    }

    pub fn process_even_order_cascade(&mut self, x: f32x8) -> f32x8 {
        self.filter2
            .process_even_order_cascade(self.filter1.process_even_order_cascade(x))
    }

    pub fn process_odd_order_cascade(&mut self, x: f32x8) -> f32x8 {
        self.filter2
            .process_odd_order_cascade(self.filter1.process_odd_order_cascade(x))
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, f32x8) -> f32x8 {
        match process_type {
            ProcessType::ProcessIIR1Only => WideF32LinkwitzRileyBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => WideF32LinkwitzRileyBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => {
                WideF32LinkwitzRileyBand::process_even_order_cascade
            }
            ProcessType::ProcessOddOrderCascade => {
                WideF32LinkwitzRileyBand::process_odd_order_cascade
            }
        }
    }

    pub fn update(&mut self, lw_coeffs: &WideF32LinkwitzRileyCoefficients) {
        self.filter1.update(&lw_coeffs.coeffs);
        self.filter2.update(&lw_coeffs.coeffs);
        self.process = WideF32LinkwitzRileyBand::get_process(lw_coeffs.coeffs.process);
    }
}
