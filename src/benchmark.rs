#[cfg(test)]
mod tests {
    #![feature(test)]
    extern crate test; //rust-analyser complains, but this should work in nightly

    use crate::filter_band::{FilterBand, FilterBandCoefficients};
    use test::Bencher;

    fn rand(x: f32) -> f32 {
        ((x * 12.9898).sin() * 43758.5453).fract()
    }

    fn static_filter_benchmark_1(
        input_samples: &Vec<f32>,
        output_samples: &mut Vec<f32>,
        order: f32,
    ) {
        let fs = 48000.0;

        let mut filter1 =
            FilterBand::from(&FilterBandCoefficients::highpass(100.0, 1.0, order, fs));
        let mut filter2 =
            FilterBand::from(&FilterBandCoefficients::lowpass(5000.0, 1.0, order, fs));
        let mut filter3 = FilterBand::from(&FilterBandCoefficients::highshelf(
            2000.0, 6.0, 1.0, order, fs,
        ));
        let mut filter4 = FilterBand::from(&FilterBandCoefficients::bell(3000.0, -6.0, 1.0, fs));

        for (input, output) in input_samples.iter().zip(output_samples.iter_mut()) {
            let mut sample = *input;
            sample = (filter1.process)(&mut filter1, sample);
            sample = (filter2.process)(&mut filter2, sample);
            sample = (filter3.process)(&mut filter3, sample);
            sample = (filter4.process)(&mut filter4, sample);
            *output = sample;
        }
    }

    #[bench]
    fn test_static_filter_benchmark_1(b: &mut Bencher) {
        let input_samples: Vec<f32> = (0..100000).map(|x| rand(x as f32)).collect();
        let mut output_samples: Vec<f32> = (0..100000).map(|_| 0.0).collect();
        b.iter(|| {
            test::black_box(static_filter_benchmark_1(
                &input_samples,
                &mut output_samples,
                3.0,
            ));
        });
    }
}
