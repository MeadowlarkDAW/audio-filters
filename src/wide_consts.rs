use wide::f32x8;
use wide::f64x4;

#[allow(non_camel_case_types)]
#[repr(C, align(16))]
union ConstUnionHack128bit {
    f64a4: [f64; 4],
    f32a8: [f32; 8],
    f32a4: [f32; 4],
    i32a4: [i32; 4],
    i64a4: [i64; 4],
    f64a2: [f64; 2],
    f32x8: f32x8,
    f64x4: f64x4,
    u128: u128,
}

macro_rules! const_f64_as_f64x4 {
    ($i:ident, $f:expr) => {
        pub const $i: f64x4 = unsafe {
            ConstUnionHack128bit {
                f64a4: [$f, $f, $f, $f],
            }
            .f64x4
        };
    };
}

macro_rules! const_f64_as_f32x8 {
    ($i:ident, $f:expr) => {
        pub const $i: f32x8 = unsafe {
            ConstUnionHack128bit {
                f64a4: [$f, $f, $f, $f],
            }
            .f32x8
        };
    };
}

#[allow(dead_code)]
#[allow(non_snake_case)]
pub(crate) mod F32x8 {
    use super::*;
    const_f64_as_f32x8!(N1, 1.0);
    const_f64_as_f32x8!(N2, 2.0);
    const_f64_as_f32x8!(N3, 3.0);
    const_f64_as_f32x8!(N4, 4.0);
    const_f64_as_f32x8!(N5, 5.0);
    const_f64_as_f32x8!(N6, 6.0);
    const_f64_as_f32x8!(N7, 7.0);
    const_f64_as_f32x8!(N8, 8.0);
    const_f64_as_f32x8!(N9, 9.0);
    const_f64_as_f32x8!(N10, 10.0);
    const_f64_as_f32x8!(N20, 20.0);
    const_f64_as_f32x8!(N40, 40.0);
}
#[allow(dead_code)]
#[allow(non_snake_case)]
pub(crate) mod F64x4 {
    use super::*;
    const_f64_as_f64x4!(N1, 1.0);
    const_f64_as_f64x4!(N2, 2.0);
    const_f64_as_f64x4!(N3, 3.0);
    const_f64_as_f64x4!(N4, 4.0);
    const_f64_as_f64x4!(N5, 5.0);
    const_f64_as_f64x4!(N6, 6.0);
    const_f64_as_f64x4!(N7, 7.0);
    const_f64_as_f64x4!(N8, 8.0);
    const_f64_as_f64x4!(N9, 9.0);
    const_f64_as_f64x4!(N10, 10.0);
    const_f64_as_f64x4!(N20, 20.0);
    const_f64_as_f64x4!(N40, 40.0);
}
