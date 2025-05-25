use crate::algebra_error::AlgebraResult;
use crate::primitives::point::Point;
use crate::surfaces::nurbs_surface::NurbsSurface;

pub fn surface_surface_intersection_non_overlap(
    _surface1: &NurbsSurface,
    _surface2: &NurbsSurface,
) -> AlgebraResult<Vec<Point>> {
    Ok(vec![])
}


#[cfg(test)]
mod tests {
    use crate::{face::Face, primitives::{color::Color10, efloat::EFloat64, primitive_scene::PrimitiveScene}};

    use super::*;

    #[test]
    fn test_surface_surface_intersection_non_overlap() -> AlgebraResult<()> {
        let surface1 = NurbsSurface::try_new(
            vec![
                vec![
                    Point::from_f64(0.0, 2.0, 0.0),
                    Point::from_f64(1.0, 1.0, 0.0),
                    Point::from_f64(2.0, 0.0, 0.0),
                    Point::from_f64(3.0, 1.0, 0.0),
                ],
                vec![
                    Point::from_f64(0.0, 1.0, 1.0),
                    Point::from_f64(1.0, 0.0, 1.0),
                    Point::from_f64(2.0, 1.0, 1.0),
                    Point::from_f64(3.0, 0.0, 1.0),
                ],
                vec![
                    Point::from_f64(0.0, 0.0, 2.0),
                    Point::from_f64(1.0, 1.0, 2.0),
                    Point::from_f64(2.0, 0.0, 2.0),
                    Point::from_f64(3.0, 1.0, 2.0),
                ],
                vec![
                    Point::from_f64(0.0, 1.0, 3.0),
                    Point::from_f64(1.0, 0.0, 3.0),
                    Point::from_f64(2.0, 1.0, 3.0),
                    Point::from_f64(3.0, 0.0, 3.0),
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
        let face2 = Face::try_new_from_surface(surface2)?;

        let mut scene = PrimitiveScene::new();
        scene.add_face(&face1, Color10::Blue)?;
        scene.add_face(&face2, Color10::Red)?;
        scene.add_convex_hull(surface1.get_convex_hull()?, Color10::Blue);

        scene.save_to_file("test_outputs/nurbs_surface_intersection_non_overlap.html")?;

        Ok(())
    }
}

