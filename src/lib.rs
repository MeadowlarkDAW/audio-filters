#![no_std]

pub mod consts;
pub mod filter_band;
mod first_order_iir;
pub mod second_order_iir;
pub mod units;

const MAX_POLE_COUNT: usize = 64;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
