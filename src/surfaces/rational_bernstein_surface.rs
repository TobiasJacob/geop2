use crate::{
    bernstein::bernstein_surface::BernsteinSurface,
    primitives::{efloat::EFloat64, point::Point},
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
