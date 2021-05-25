# Audio Filters - WIP
A collection of filters for real-time audio processing

The overall design is in progress, there can be frequent (usually minor) breaking changes.

## Feature Progress
    
- [x] `#![no_std]` (via [libm](https://github.com/rust-lang/libm))
- [x] f32 & f64 capable (via [num-traits](https://github.com/rust-num/num-traits))
- [ ] SIMD (Some experimental support)
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
- [x] Linkwitz-Riley filters
- [ ] Elliptic filters
- [ ] Phase aligned crossovers
- [ ] Decramping near nyquist
- [x] Minimum Phase IIR Mode
- [ ] Linear Phase Mode

```rust
let sample_rate_hz = 48000.0;
let cutoff_hz = 1000.0;
let gain_db = 6.0;
let width_oct = 1.0;
let order = 4.0;

let coeffs = FilterBandCoefficients::highshelf(cutoff_hz, gain_db, width_oct, order, sample_rate_hz);

let mut filter_left = FilterBand::from(&coeffs);
let mut filter_right = FilterBand::from(&coeffs);

for i in 0..1000 {
    left[i] = (filter_left.process)(&mut filter_left, left[i]);
    right[i] = (filter_right.process)(&mut filter_right, right[i]);
}
```

