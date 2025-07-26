use std::fmt::Display;

use crate::{
    algebra_error::AlgebraResult,
    primitives::{convex_hull::ConvexHull, efloat::EFloat64, point::Point},
};

pub trait SurfaceLike: Display {
    /// Returns the valid parameter span for u (min_u, max_u).
    fn u_span(&self) -> (EFloat64, EFloat64);

    /// Returns the valid parameter span for v (min_v, max_v).
    fn v_span(&self) -> (EFloat64, EFloat64);

    /// Evaluates the surface at the parameter pair (u, v).
    fn eval(&self, u: EFloat64, v: EFloat64) -> Point;

    /// Subdivides the surface along the u–direction at parameter `t` into two new surface segments.
    fn subdivide_u(&self, t: EFloat64) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized;

    /// Subdivides the surface along the v–direction at parameter `t` into two new surface segments.
    fn subdivide_v(&self, t: EFloat64) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized;

    fn split_u(&self) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized {
        let u_range = self.u_span().1 - self.u_span().0;
        let u = self.u_span().0 + (u_range / EFloat64::from(2.0)).unwrap_or(EFloat64::zero());
        self.subdivide_u(u)
    }
    fn split_v(&self) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized {
        let v_range = self.v_span().1 - self.v_span().0;
        let v = self.v_span().0 + (v_range / EFloat64::from(2.0)).unwrap_or(EFloat64::zero());
        self.subdivide_v(v)
    }

    /// Returns the convex hull of the control points of the surface.
    /// This is the smallest convex set containing all control points.
    fn get_convex_hull(&self) -> AlgebraResult<ConvexHull>;
}
