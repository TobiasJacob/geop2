use crate::{
    primitives::{efloat::EFloat64, point::Point},
    surfaces::bernstein_surface::BernsteinSurface,
};

pub struct RationalBernsteinSurface {
    pub nominator: BernsteinSurface<Point>,
    pub denominator: BernsteinSurface<EFloat64>,
}

impl RationalBernsteinSurface {
    pub fn new(
        nominator: BernsteinSurface<Point>,
        denominator: BernsteinSurface<EFloat64>,
    ) -> Self {
        Self {
            nominator,
            denominator,
        }
    }
}
