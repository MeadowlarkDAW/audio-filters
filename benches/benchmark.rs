#![feature(test)]
extern crate test; //rust-analyser complains, but this should work in nightly

use test::Bencher;
use wide::{f32x8, f64x4};
use {
    audio_filters::filter_band::{FilterBand, FilterBandCoefficients},
    audio_filters::filter_band_wide::{
        WideF32FilterBand, WideF32FilterBandCoefficients, WideF64FilterBand,
        WideF64FilterBandCoefficients,
    },
};

fn rand(x: f32) -> f32 {
    ((x * 12.9898).sin() * 43758.5453).fract()
}

fn rand64(x: f64) -> f64 {
    ((x * 12.98986578567).sin() * 43758.54535678567).fract()
}

fn dynamic_filter_benchmark_1(input_samples: &Vec<f32>, output_samples: &mut Vec<f32>, order: f32) {
    let fs = 48000.0;

    let mut filter1 = FilterBand::from(&FilterBandCoefficients::highpass(100.0, 1.0, order, fs));
    let mut filter2 = FilterBand::from(&FilterBandCoefficients::lowpass(5000.0, 1.0, order, fs));
    let mut filter3 = FilterBand::from(&FilterBandCoefficients::highshelf(
        2000.0, 6.0, 1.0, order, fs,
    ));
    let mut filter4 = FilterBand::from(&FilterBandCoefficients::bell(3000.0, -6.0, 1.0, fs));

    for (i, (input, output)) in input_samples
        .iter()
        .zip(output_samples.iter_mut())
        .enumerate()
    {
        let mut sample = *input;
        sample = (filter1.process)(&mut filter1, sample);
        sample = (filter2.process)(&mut filter2, sample);
        sample = (filter3.process)(&mut filter3, sample);
        sample = (filter4.process)(&mut filter4, sample);
        *output = sample;
        let n = (i as f32).sin();
        filter1.update(&FilterBandCoefficients::highpass(
            n * 100.0 + 100.0,
            1.0,
            order,
            fs,
        ));
        filter2.update(&FilterBandCoefficients::lowpass(
            n * 100.0 + 2000.0,
            1.0,
            order,
            fs,
        ));
        filter3.update(&FilterBandCoefficients::highshelf(
            n * 100.0 + 1000.0,
            n,
            1.0,
            order,
            fs,
        ));
        filter4.update(&FilterBandCoefficients::bell(n * 100.0 + 500.0, n, 1.0, fs));
    }
}

fn static_filter_benchmark_1(input_samples: &Vec<f32>, output_samples: &mut Vec<f32>, order: f32) {
    let fs = 48000.0;

    let mut filter1 = FilterBand::from(&FilterBandCoefficients::highpass(100.0, 1.0, order, fs));
    let mut filter2 = FilterBand::from(&FilterBandCoefficients::lowpass(5000.0, 1.0, order, fs));
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
    let input_samples: Vec<f32> = (0..1000000).map(|x| rand(x as f32)).collect();
    let mut output_samples: Vec<f32> = (0..1000000).map(|_| 0.0).collect();
    b.iter(|| {
        test::black_box(static_filter_benchmark_1(
            &input_samples,
            &mut output_samples,
            9.0,
        ));
    });
}

#[bench]
fn test_dynamic_filter_benchmark_1(b: &mut Bencher) {
    let input_samples: Vec<f32> = (0..100000).map(|x| rand(x as f32)).collect();
    let mut output_samples: Vec<f32> = (0..100000).map(|_| 0.0).collect();
    b.iter(|| {
        test::black_box(dynamic_filter_benchmark_1(
            &input_samples,
            &mut output_samples,
            9.0,
        ));
    });
}

fn static_filter_benchmark_1_wide32x8(
    input_samples: &Vec<f32x8>,
    output_samples: &mut Vec<f32x8>,
    order: f32,
) {
    let fs = 48000.0;

    let mut filter1 = WideF32FilterBand::from(&WideF32FilterBandCoefficients::from(
        FilterBandCoefficients::highpass(100.0, 1.0, order, fs),
    ));
    let mut filter2 = WideF32FilterBand::from(&WideF32FilterBandCoefficients::from(
        FilterBandCoefficients::lowpass(5000.0, 1.0, order, fs),
    ));
    let mut filter3 = WideF32FilterBand::from(&WideF32FilterBandCoefficients::from(
        FilterBandCoefficients::highshelf(2000.0, 6.0, 1.0, order, fs),
    ));
    let mut filter4 = WideF32FilterBand::from(&WideF32FilterBandCoefficients::from(
        FilterBandCoefficients::bell(3000.0, -6.0, 1.0, fs),
    ));

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
fn test_static_filter_benchmark_1_wide32x8(b: &mut Bencher) {
    let inputs = (0..1000000)
        .map(|x| {
            f32x8::from([
                rand((x * 1) as f32),
                rand((x * 2) as f32),
                rand((x * 3) as f32),
                rand((x * 4) as f32),
                rand((x * 5) as f32),
                rand((x * 6) as f32),
                rand((x * 7) as f32),
                rand((x * 8) as f32),
            ])
        })
        .collect();

    let mut outputs = (0..1000000).map(|_| f32x8::ZERO).collect();

    b.iter(|| {
        test::black_box(static_filter_benchmark_1_wide32x8(
            &inputs,
            &mut outputs,
            9.0,
        ));
    });
}

fn static_filter_benchmark_1_wide64x4(
    input_samples: &Vec<f64x4>,
    output_samples: &mut Vec<f64x4>,
    order: f64,
) {
    let fs = 48000.0;

    let mut filter1 = WideF64FilterBand::from(&WideF64FilterBandCoefficients::from(
        FilterBandCoefficients::highpass(100.0, 1.0, order, fs),
    ));
    let mut filter2 = WideF64FilterBand::from(&WideF64FilterBandCoefficients::from(
        FilterBandCoefficients::lowpass(5000.0, 1.0, order, fs),
    ));
    let mut filter3 = WideF64FilterBand::from(&WideF64FilterBandCoefficients::from(
        FilterBandCoefficients::highshelf(2000.0, 6.0, 1.0, order, fs),
    ));
    let mut filter4 = WideF64FilterBand::from(&WideF64FilterBandCoefficients::from(
        FilterBandCoefficients::bell(3000.0, -6.0, 1.0, fs),
    ));

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
fn test_static_filter_benchmark_1_wide64x4(b: &mut Bencher) {
    let inputs = (0..1000000)
        .map(|x| {
            f64x4::from([
                rand64((x * 1) as f64),
                rand64((x * 2) as f64),
                rand64((x * 3) as f64),
                rand64((x * 4) as f64),
            ])
        })
        .collect();

    let mut outputs = (0..1000000).map(|_| f64x4::ZERO).collect();

    b.iter(|| {
        test::black_box(static_filter_benchmark_1_wide64x4(
            &inputs,
            &mut outputs,
            9.0,
        ));
    });
}
