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
let fs = 48000.0;
let f0 = 1000.0;
let gain = 6.0;
let width = 1.0;
let slope = 4.0;

let coeffs = FilterBandCoefficients::highshelf(f0, gain, width, slope, fs);

let mut filter_left = FilterBand::from(&coeffs);
let mut filter_right = FilterBand::from(&coeffs);

for i in 0..1000 {
    left[i] = (filter_left.process)(&mut filter_left, left[i]);
    right[i] = (filter_right.process)(&mut filter_right, right[i]);
}
```

