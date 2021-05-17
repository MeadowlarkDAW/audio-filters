# Audio Filters
A collection of filters for real-time audio processing

## Feature Progress
    
- [x] `#![no_std]` (via [libm](https://github.com/rust-lang/libm))
- [x] f32 & f64 capable (via [num-traits](https://github.com/rust-num/num-traits))
- [ ] SIMD
- [ ] Documentation
- [ ] Tests

### Filter Types

- [x] Bell
- [x] Low Pass
- [x] High Pass
- [x] Low Shelf
- [x] High Shelf
- [x] Notch
- [x] Band Pass
- [x] All Pass
- [ ] Higher Order Bell
- [ ] Higher Order Band Pass
- [ ] Asymmetrical
- [ ] Tilt
- [ ] Flat Tilt

### Filter Features

- [x] Bode Plot (phase & amplitude)
- [x] 1st and 2nd order filter primitives
- [x] Higher orders via cascading
- [x] Virtual analog (VA) State Variable Filters (SVF) for both 1st & 2nd order IIR.
- [ ] Linkwitz-Riley filters
- [ ] Elliptic filters
- [ ] Phase aligned crossovers
- [ ] Decramping near nyquist
- [x] Minimum Phase IIR Mode
- [ ] Linear Phase Mode

```rust
let sample_rate = 48000.0;
let f0 = 1000.0;
let gain = 6.0;
let bandwidth = 1.0;
let slope = 4.0;

let mut filter = StereoFilterBand::new(sample_rate);

filter.update(BandType::HighShelf, f0, gain, bandwidth, slope, sample_rate);

for i in 0..1000 {
    let [l_out, r_out] = filter.process(left[i], right[i]);
    left[i] = l_out;
    right[i] = r_out;
}
```

