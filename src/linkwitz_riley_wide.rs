use crate::{
    filter_band::ProcessType,
    filter_band_wide::{WideFilterBand, WideFilterBandCoefficients},
    wide_units::WIDE,
};

#[derive(Copy, Clone, Debug)]
pub struct WideLinkwitzRileyCoefficients<T: WIDE> {
    pub coeffs: WideFilterBandCoefficients<T>,
}

impl<T: WIDE> WideLinkwitzRileyCoefficients<T> {
    //The resulting Linkwitz-Riley filter will have 2x to order of the input coefficients and 2x gain
    pub fn from(coeffs: WideFilterBandCoefficients<T>) -> Self {
        WideLinkwitzRileyCoefficients { coeffs }
    }
}

#[derive(Copy, Clone)]
pub struct WideLinkwitzRileyBand<T: WIDE> {
    pub filter1: WideFilterBand<T>,
    pub filter2: WideFilterBand<T>,
    pub process: fn(&mut Self, T) -> T,
}

impl<T: WIDE> WideLinkwitzRileyBand<T> {
    pub fn from(lw_coeffs: &WideLinkwitzRileyCoefficients<T>) -> WideLinkwitzRileyBand<T> {
        WideLinkwitzRileyBand {
            filter1: WideFilterBand::from(&lw_coeffs.coeffs),
            filter2: WideFilterBand::from(&lw_coeffs.coeffs),
            process: WideLinkwitzRileyBand::get_process(lw_coeffs.coeffs.process),
        }
    }

    pub fn process_iir1_only(&mut self, x: T) -> T {
        self.filter2
            .process_iir1_only(self.filter1.process_iir1_only(x))
    }

    pub fn process_iir2_only(&mut self, x: T) -> T {
        self.filter2
            .process_iir2_only(self.filter1.process_iir2_only(x))
    }

    pub fn process_even_order_cascade(&mut self, x: T) -> T {
        self.filter2
            .process_even_order_cascade(self.filter1.process_even_order_cascade(x))
    }

    pub fn process_odd_order_cascade(&mut self, x: T) -> T {
        self.filter2
            .process_odd_order_cascade(self.filter1.process_odd_order_cascade(x))
    }

    pub fn get_process(process_type: ProcessType) -> fn(&mut Self, T) -> T {
        match process_type {
            ProcessType::ProcessIIR1Only => WideLinkwitzRileyBand::process_iir1_only,
            ProcessType::ProcessIIR2Only => WideLinkwitzRileyBand::process_iir2_only,
            ProcessType::ProcessEvenOrderCascade => {
                WideLinkwitzRileyBand::process_even_order_cascade
            }
            ProcessType::ProcessOddOrderCascade => WideLinkwitzRileyBand::process_odd_order_cascade,
        }
    }

    pub fn update(&mut self, lw_coeffs: &WideLinkwitzRileyCoefficients<T>) {
        self.filter1.update(&lw_coeffs.coeffs);
        self.filter2.update(&lw_coeffs.coeffs);
        self.process = WideLinkwitzRileyBand::get_process(lw_coeffs.coeffs.process);
    }
}
