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
                "Creating traingle with points: {:?}, {:?}, {:?}",
                &a, &b, &c
            ))
        };

        if a == b || b == c || c == a {
            return Err(AlgebraError::new("a, b, c must be distinct points".into()))
                .with_context(&context);
        }

        // Calculate the normal direction using the cross product.
        let normal = (b - a).cross(c - a);
        // let normal2 = (c - b).cross(a - b);
        // let normal3 = (a - c).cross(b - c);
        // let mut normal = normal1;
        // if normal2.norm_sq().lower_bound > normal.norm_sq().lower_bound {
        //     normal = normal2;
        // }
        // if normal3.norm_sq().lower_bound > normal.norm_sq().lower_bound {
        //     normal = normal3;
        // }
        // if normal.is_zero() {
        //     return Err(AlgebraError::new(
        //         format!(
        //             "Cannot find a normal for this triangle. Candidates are {:?}, {:?}, {:?}",
        //             normal1, normal2, normal3
        //         )
        //         .into(),
        //     ))
        //     .with_context(&context);
        // }

        let normal = normal.normalize().with_context(&context)?;
        Ok(TriangleFace { a, b, c, normal })
    }

    pub fn try_new_with_normal(
        a: Point,
        b: Point,
        c: Point,
        normal: Point,
    ) -> AlgebraResult<TriangleFace> {
        let context = |err: AlgebraError| {
            err.with_context(format!(
                "Creating traingle with points: {}, {}, {}",
                &a, &b, &c
            ))
        };

        if normal.is_zero() {
            return Err(AlgebraError::new("normal must be a non-zero vector".into()))
                .with_context(&context);
        }

        // check that normal is perpendicular to ab and ac and bc
        let ab = b - a;
        let ac = c - a;
        let bc = c - b;
        if normal.dot(ab) != 0.0 || normal.dot(ac) != 0.0 || normal.dot(bc) != 0.0 {
            return Err(AlgebraError::new(
                "normal must be perpendicular to ab and ac and bc".into(),
            ))
            .with_context(&context);
        }

        // check normal is normalized
        if normal.norm() != 1.0 {
            return Err(AlgebraError::new("normal must be normalized".into()))
                .with_context(&context);
        }

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
