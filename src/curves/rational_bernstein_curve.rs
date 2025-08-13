use crate::primitives::{efloat::EFloat64, point::Point};

pub struct RationalBernsteinCurve {
    pub control_points: Vec<Point>,
    pub weights: Vec<EFloat64>,
}

impl RationalBernsteinCurve {
    pub fn new(control_points: Vec<Point>, weights: Vec<EFloat64>) -> Self {
        Self {
            control_points,
            weights,
        }
    }
}
