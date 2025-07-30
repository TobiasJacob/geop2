use crate::algebra_error::AlgebraResult;
use crate::face::Face;
use crate::primitives::efloat::EFloat64;
use crate::primitives::point::Point;
use crate::surfaces::nurbs_surface::NurbsSurface;
use crate::surfaces::surface_like::SurfaceLike;

// After this, refinement will be done, but no more points will be added to the curve.
const INTERSECTION_SMALLEST_PATCH: f64 = 1e-2;

// Tries to find a point that lies on both surfaces and is as accurate as possible.
// If there are multiple such points, it returns an arbitrary one.
// In the special case that it intersects with the boundary of the surface, it returns the point on the boundary.
// If there is a point with colinear normals, it returns this point.
fn _surface_surface_intersection_non_overlap_refinement(
    _surface1: &NurbsSurface,
    _surface2: &NurbsSurface,
) -> AlgebraResult<Point> {
    Ok(Point::zero())
}

// 
pub fn surface_surface_intersection_non_overlap(
    surface1: &NurbsSurface,
    surface2: &NurbsSurface,
    divide_u_or_v: bool
) -> AlgebraResult<Vec<NurbsSurface>> {
    if surface1
        .get_convex_hull()?
        .intersects(&surface2.get_convex_hull()?) == false
    {
        return Ok(vec![]);
    }

    if surface1.u_span().1 - surface1.u_span().0 <= EFloat64::from(INTERSECTION_SMALLEST_PATCH) {
        return Ok(vec![surface1.clone(), surface2.clone()]);
    }

    let (s1a, s1b) = if divide_u_or_v {
        surface1.split_u()?
    } else {
        surface1.split_v()?
    };
    let (s2a, s2b) = if divide_u_or_v {
        surface2.split_u()?
    } else {
        surface2.split_v()?
    };

    let intersections = vec![
        surface_surface_intersection_non_overlap(&s1a, &s2a, !divide_u_or_v)?,
        surface_surface_intersection_non_overlap(&s1a, &s2b, !divide_u_or_v)?,
        surface_surface_intersection_non_overlap(&s1b, &s2a, !divide_u_or_v)?,
        surface_surface_intersection_non_overlap(&s1b, &s2b, !divide_u_or_v)?,
    ].iter().flatten()
        .map(|s| s.clone())
        .collect::<Vec<NurbsSurface>>();

    return Ok(intersections);
}

pub fn face_face_intersection_special_points(
    _face1: &Face,
    _face2: &Face,
) -> AlgebraResult<Vec<Point>> {
    todo!("Implement intersection for special points in face_face_intersection_special_points")
    // let mut result = vec![];
    // for b in face1.bounds {
    //     for c in b.bounds {
    //     }
    // }
    // Ok(result)
}


#[cfg(test)]
mod tests {
    use crate::{face::Face, intersection::test_surface_1, primitives::{color::Color10, efloat::EFloat64, primitive_scene::PrimitiveScene}};

    use super::*;

    #[test]
    fn test_surface_surface_intersection_non_overlap() -> AlgebraResult<()> {
        let surface1 = test_surface_1()?;
        let face1 = Face::try_new_from_surface(surface1.clone())?;

        let surface2 = NurbsSurface::try_new(
            vec![
                vec![
                    Point::from_f64(0.0, 1.0, 0.0),
                    Point::from_f64(1.0, 4.0, 0.0),
                    Point::from_f64(2.0, 3.0, 0.0),
                    Point::from_f64(3.0, 0.0, 0.0),
                ],
                vec![
                    Point::from_f64(0.0, 0.0, 1.0),
                    Point::from_f64(1.0, 1.0, 1.0),
                    Point::from_f64(2.0, 0.0, 1.0),
                    Point::from_f64(3.0, 1.0, 1.0),
                ],
                vec![
                    Point::from_f64(0.0, 1.0, 2.0),
                    Point::from_f64(1.0, 0.0, 2.0),
                    Point::from_f64(2.0, 1.0, 2.0),
                    Point::from_f64(3.0, 0.0, 2.0),
                ],
                vec![
                    Point::from_f64(0.0, 0.0, 3.0),
                    Point::from_f64(1.0, 1.0, 3.0),
                    Point::from_f64(2.0, 0.0, 3.0),
                    Point::from_f64(3.0, 1.0, 3.0),
                ],
            ],
            vec![
                vec![
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                ],
                vec![
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                ],
                vec![
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                ],
                vec![
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                ],
            ],
            vec![
                EFloat64::from(0.0),
                EFloat64::from(0.0),
                EFloat64::from(0.0),
                EFloat64::from(0.5),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
            ],
            vec![
                EFloat64::from(0.0),
                EFloat64::from(0.0),
                EFloat64::from(0.0),
                EFloat64::from(0.5),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
            ],
            2,
            2,
        )?;
        let face2 = Face::try_new_from_surface(surface2.clone())?;

        let mut scene = PrimitiveScene::new();
        scene.add_face_wireframe(&face1, Color10::Blue, 20)?;
        scene.add_face_wireframe(&face2, Color10::Red, 20)?;

        let intersection = surface_surface_intersection_non_overlap(&surface1, &surface2, true)?;
        for surface in intersection {
            let face = Face::try_new_from_surface(surface)?;
            scene.add_face(&face, Color10::Green, 1)?;
        }

        scene.save_to_file("test_outputs/nurbs_surface_intersection_non_overlap.html")?;

        Ok(())
    }
}

