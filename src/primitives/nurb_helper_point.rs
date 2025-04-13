use crate::algebra_error::AlgebraResult;

use super::{efloat::EFloat64, point::Point};

/// Helper struct for representing homogeneous control points during evaluation.
#[derive(Clone)]
pub struct NurbHelperPoint {
    point: Point,
    weight: EFloat64,
}

impl NurbHelperPoint {
    pub fn new(point: Point, weight: EFloat64) -> Self {
        Self { point, weight }
    }

    pub fn zero() -> Self {
        Self {
            point: Point::zero(),
            weight: EFloat64::zero(),
        }
    }

    pub fn to_point(&self) -> AlgebraResult<Point> {
        self.point / self.weight
    }
}

impl std::ops::Add for NurbHelperPoint {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            point: self.point + other.point,
            weight: self.weight + other.weight,
        }
    }
}

impl std::ops::Mul<EFloat64> for NurbHelperPoint {
    type Output = Self;
    fn mul(self, scalar: EFloat64) -> Self {
        Self {
            point: self.point * scalar,
            weight: self.weight * scalar,
        }
    }
}

impl std::ops::Mul<NurbHelperPoint> for EFloat64 {
    type Output = NurbHelperPoint;
    fn mul(self, rhs: NurbHelperPoint) -> NurbHelperPoint {
        rhs * self
    }
}
