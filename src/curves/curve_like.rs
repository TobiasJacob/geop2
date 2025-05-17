use crate::algebra_error::AlgebraResult;
use crate::primitives::convex_hull::ConvexHull;
use crate::primitives::efloat::EFloat64;
use crate::primitives::point::Point;

pub trait CurveLike {
    /// Evaluates the curve at parameter t.
    fn eval(&self, t: EFloat64) -> Point;

    /// Subdivides the curve at parameter t into two new curve segments.
    fn subdivide(&self, t: EFloat64) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized;

    /// Generates the convex hull of the control polygon.
    fn control_polygon_hull(&self) -> AlgebraResult<ConvexHull>;
}
