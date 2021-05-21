use core::ops::{Add, Div, Mul, Sub};

use num_traits::NumCast;
use wide::f32x4;
use wide::f32x8;
use wide::f64x2;
use wide::f64x4;

use crate::units::FP;

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
    f32x4: f32x4,
    f64x2: f64x2,
    f64x4: f64x4,
    u128: u128,
}

macro_rules! const_f64_as_f64x4 {
    ($i:ident, $f:expr) => {
        const $i: f64x4 = unsafe {
            ConstUnionHack128bit {
                f64a4: [$f, $f, $f, $f],
            }
            .f64x4
        };
    };
}

macro_rules! const_f64_as_f64x2 {
    ($i:ident, $f:expr) => {
        const $i: f64x2 = unsafe { ConstUnionHack128bit { f64a2: [$f, $f] }.f64x2 };
    };
}

macro_rules! const_f32_as_f32x8 {
    ($i:ident, $f:expr) => {
        const $i: f32x8 = unsafe {
            ConstUnionHack128bit {
                f32a8: [$f, $f, $f, $f, $f, $f, $f, $f],
            }
            .f32x8
        };
    };
}

macro_rules! const_f32_as_f32x4 {
    ($i:ident, $f:expr) => {
        const $i: f32x4 = unsafe {
            ConstUnionHack128bit {
                f32a4: [$f, $f, $f, $f],
            }
            .f32x4
        };
    };
}

pub trait WIDE:
    Sized
    + Copy
    + Add
    + Mul
    + Sub
    + Mul<Self, Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Div<Self, Output = Self>
{
    fn from_w<T: FP>(n: T) -> Self;
    const ZERO: Self;
    const ONE: Self;
    const N0: Self;
    const N00_5: Self;
    const N0_5: Self;
    const N1: Self;
    const N2: Self;
    const N3: Self;
    const N4: Self;
    const N5: Self;
    const N6: Self;
    const N7: Self;
    const N8: Self;
    const N9: Self;
    const N10: Self;
    const N20: Self;
    const N40: Self;
}

impl WIDE for f64x4 {
    #[inline]
    fn from_w<T: FP>(n: T) -> f64x4 {
        let n: f64 = NumCast::from(n).unwrap();
        Self::from(n)
    }
    const ZERO: f64x4 = f64x4::ZERO;
    const ONE: f64x4 = f64x4::ONE;

    const_f64_as_f64x4!(N0, 0.0);
    const_f64_as_f64x4!(N00_5, 0.05);
    const_f64_as_f64x4!(N0_5, 0.5);
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
impl WIDE for f64x2 {
    #[inline]
    fn from_w<T: FP>(n: T) -> f64x2 {
        let n: f64 = NumCast::from(n).unwrap();
        Self::from(n)
    }
    const ZERO: f64x2 = f64x2::ZERO;
    const ONE: f64x2 = f64x2::ONE;

    const_f64_as_f64x2!(N0, 0.0);
    const_f64_as_f64x2!(N00_5, 0.05);
    const_f64_as_f64x2!(N0_5, 0.5);
    const_f64_as_f64x2!(N1, 1.0);
    const_f64_as_f64x2!(N2, 2.0);
    const_f64_as_f64x2!(N3, 3.0);
    const_f64_as_f64x2!(N4, 4.0);
    const_f64_as_f64x2!(N5, 5.0);
    const_f64_as_f64x2!(N6, 6.0);
    const_f64_as_f64x2!(N7, 7.0);
    const_f64_as_f64x2!(N8, 8.0);
    const_f64_as_f64x2!(N9, 9.0);
    const_f64_as_f64x2!(N10, 10.0);
    const_f64_as_f64x2!(N20, 20.0);
    const_f64_as_f64x2!(N40, 40.0);
}
impl WIDE for f32x8 {
    #[inline]
    fn from_w<T: FP>(n: T) -> f32x8 {
        let n: f32 = NumCast::from(n).unwrap();
        Self::from(n)
    }
    const ZERO: f32x8 = f32x8::ZERO;
    const ONE: f32x8 = f32x8::ONE;
    const_f32_as_f32x8!(N0, 0.0);
    const_f32_as_f32x8!(N00_5, 0.05);
    const_f32_as_f32x8!(N0_5, 0.5);
    const_f32_as_f32x8!(N1, 1.0);
    const_f32_as_f32x8!(N2, 2.0);
    const_f32_as_f32x8!(N3, 3.0);
    const_f32_as_f32x8!(N4, 4.0);
    const_f32_as_f32x8!(N5, 5.0);
    const_f32_as_f32x8!(N6, 6.0);
    const_f32_as_f32x8!(N7, 7.0);
    const_f32_as_f32x8!(N8, 8.0);
    const_f32_as_f32x8!(N9, 9.0);
    const_f32_as_f32x8!(N10, 10.0);
    const_f32_as_f32x8!(N20, 20.0);
    const_f32_as_f32x8!(N40, 40.0);
}

impl WIDE for f32x4 {
    #[inline]
    fn from_w<T: FP>(n: T) -> f32x4 {
        let n: f32 = NumCast::from(n).unwrap();
        Self::from(n)
    }
    const ZERO: f32x4 = f32x4::ZERO;
    const ONE: f32x4 = f32x4::ONE;
    const_f32_as_f32x4!(N0, 0.0);
    const_f32_as_f32x4!(N00_5, 0.05);
    const_f32_as_f32x4!(N0_5, 0.5);
    const_f32_as_f32x4!(N1, 1.0);
    const_f32_as_f32x4!(N2, 2.0);
    const_f32_as_f32x4!(N3, 3.0);
    const_f32_as_f32x4!(N4, 4.0);
    const_f32_as_f32x4!(N5, 5.0);
    const_f32_as_f32x4!(N6, 6.0);
    const_f32_as_f32x4!(N7, 7.0);
    const_f32_as_f32x4!(N8, 8.0);
    const_f32_as_f32x4!(N9, 9.0);
    const_f32_as_f32x4!(N10, 10.0);
    const_f32_as_f32x4!(N20, 20.0);
    const_f32_as_f32x4!(N40, 40.0);
}
