use crate::primitives::{efloat::EFloat64, point::Point};

pub trait Zero {
    fn zero() -> Self;
}

impl Zero for EFloat64 {
    fn zero() -> Self {
        EFloat64::zero()
    }
}

impl Zero for Point {
    fn zero() -> Self {
        Point::zero()
    }
}
