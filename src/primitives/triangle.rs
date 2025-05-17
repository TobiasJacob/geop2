use std::fmt::Display;

use crate::algebra_error::{AlgebraError, WithContext};
use crate::primitives::point::Point;
use crate::{algebra_error::AlgebraResult, primitives::efloat::EFloat64};

#[derive(Debug, Clone)]
pub struct TriangleFace {
    pub a: Point,
    pub b: Point,
    pub c: Point,
    pub normal: Point,
}

impl TriangleFace {
    pub fn try_new(a: Point, b: Point, c: Point) -> AlgebraResult<TriangleFace> {
        let context = |err: AlgebraError| {
            err.with_context(format!(
                "Creating traingle with points: {}, {}, {}",
                &a, &b, &c
            ))
        };

        if a == b || b == c || c == a {
            return Err(AlgebraError::new("a, b, c must be distinct points".into()))
                .with_context(&context);
        }

        // Calculate the normal direction using the cross product.
        let ab = b - a;
        let ac = c - a;
        let normal = ab.cross(ac).normalize().with_context(&context)?;
        Ok(TriangleFace { a, b, c, normal })
    }

    /// Computes the oriented distance of a point from the plane defined by this face.
    pub fn distance_to_point(&self, p: &Point) -> EFloat64 {
        let ap = *p - self.a;
        ap.dot(self.normal)
    }
}

impl Display for TriangleFace {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Face: a={}, b={}, c={}, normal={}",
            self.a, self.b, self.c, self.normal
        )
    }
}
