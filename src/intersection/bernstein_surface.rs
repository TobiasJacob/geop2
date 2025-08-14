use crate::{
    algebra_error::AlgebraResult,
    primitives::point::Point,
    surfaces::{bernstein_surface::BernsteinSurface, surface_like::SurfaceLike},
};

pub fn get_conormal_points(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
) -> AlgebraResult<Vec<Point>> {
    let hull_s1 = s1.get_convex_hull()?;
    let hull_s2 = s2.get_convex_hull()?;

    if !hull_s1.intersects(&hull_s2) {
        return Ok(vec![]);
    }

    // Now get normals
    let n1 = s1.derivative_u().cross(&s1.derivative_v());
    let n2 = s2.derivative_u().cross(&s2.derivative_v());

    let _hull_n1 = n1.get_convex_hull()?;
    let _hull_n2 = n2.get_convex_hull()?;

    // Check if the dot product

    Ok(vec![])
}
