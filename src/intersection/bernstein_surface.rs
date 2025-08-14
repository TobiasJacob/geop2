use crate::{
    algebra_error::AlgebraResult,
    bernstein::{bernstein_hypervolume::BernsteinHyperVolume, bernstein_surface::BernsteinSurface},
    primitives::{efloat::EFloat64, point::Point},
    surfaces::surface_like::SurfaceLike,
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
    let n1 = s1.normal();
    let n2 = s2.normal();

    let n1_hyper = n1.to_hypervolume_tu();
    let n2_hyper = n2.to_hypervolume_vw();

    let cross = n1_hyper.cross(&n2_hyper);

    // Check if the dot product
    let _hull = cross.get_convex_hull()?;

    Ok(vec![])
}

fn _refine(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
    resultant: &BernsteinHyperVolume<EFloat64>,
) -> AlgebraResult<Vec<Point>> {
    if !s1.get_convex_hull()?.intersects(&s2.get_convex_hull()?) {
        return Ok(vec![]);
    }

    if !resultant.could_be_zero() {
        return Ok(vec![]);
    }

    todo!("Subdivision step")
}

pub fn get_critical_points_perpendicular_tangent_x(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
) -> AlgebraResult<Vec<Point>> {
    let hull_s1 = s1.get_convex_hull()?;
    let hull_s2 = s2.get_convex_hull()?;

    if !hull_s1.intersects(&hull_s2) {
        return Ok(vec![]);
    }

    // Now get normals
    let n1 = s1.normal();
    let n2 = s2.normal();

    let n1_hyper = n1.to_hypervolume_tu();
    let n2_hyper = n2.to_hypervolume_vw();

    let intersection_tangent = n1_hyper.cross(&n2_hyper);
    let _product = intersection_tangent.dot(&BernsteinHyperVolume::unit_x());

    Ok(vec![])
}

#[cfg(test)]
mod tests {
    use crate::primitives::{color::Color10, primitive_scene::PrimitiveScene};

    use super::*;

    fn mountain() -> BernsteinSurface<Point> {
        BernsteinSurface::<Point>::new(vec![
            vec![
                Point::from_f64(0.0, 0.0, 0.0),
                Point::from_f64(1.0, 0.0, 0.0),
                Point::from_f64(2.0, 0.0, 0.0),
            ],
            vec![
                Point::from_f64(0.0, 1.0, 0.0),
                Point::from_f64(1.0, 1.0, 5.0),
                Point::from_f64(2.0, 1.0, 0.0),
            ],
            vec![
                Point::from_f64(0.0, 2.0, 0.0),
                Point::from_f64(1.0, 2.0, 0.0),
                Point::from_f64(2.0, 2.0, 0.0),
            ],
        ])
    }

    fn upside_down_mountain() -> BernsteinSurface<Point> {
        BernsteinSurface::<Point>::new(vec![
            vec![
                Point::from_f64(0.0, 0.0, 2.0),
                Point::from_f64(1.0, 0.0, 2.0),
                Point::from_f64(2.0, 0.0, 2.0),
            ],
            vec![
                Point::from_f64(0.0, 1.0, 2.0),
                Point::from_f64(1.0, 1.0, -5.0),
                Point::from_f64(2.0, 1.0, 2.0),
            ],
            vec![
                Point::from_f64(0.0, 2.0, 2.0),
                Point::from_f64(1.0, 2.0, 2.0),
                Point::from_f64(2.0, 2.0, 2.0),
            ],
        ])
    }

    #[test]
    fn test_get_distance_hypervolume() -> AlgebraResult<()> {
        let s1 = mountain();
        let s2 = upside_down_mountain();

        let distance = s1.to_hypervolume_tu() - s2.to_hypervolume_vw();
        let hull = distance.get_convex_hull().unwrap();

        let mut scene = PrimitiveScene::new();
        scene.add_surface_like(&s1, Color10::Red, 30)?;
        scene.add_surface_like(&s2, Color10::Blue, 30)?;
        scene.add_convex_hull(hull, Color10::Green);

        for c1 in distance.coefficients.iter() {
            for c2 in c1.iter() {
                for c3 in c2.iter() {
                    for c4 in c3.iter() {
                        scene.add_point(*c4, Color10::Green);
                    }
                }
            }
        }

        scene.save_to_file("test_outputs/distance_hypervolume.html")?;

        Ok(())
    }
}
