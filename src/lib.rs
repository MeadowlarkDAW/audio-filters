#![cfg_attr(not(test), no_std)]
#![feature(test)]

pub mod benchmark;
pub mod filter_band;
pub mod first_order_iir;
pub mod second_order_iir;
pub mod units;

const MAX_POLE_COUNT: usize = 32;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
