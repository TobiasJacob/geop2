use crate::{
    algebra_error::{AlgebraError, AlgebraResult},
    bernstein::{bernstein_hypervolume::BernsteinHyperVolume, bernstein_surface::BernsteinSurface},
    primitives::{
        color::Color10, efloat::EFloat64, point::Point, primitive_scene::PrimitiveScene,
        primitive_scene_recorder::PrimitiveSceneRecorder,
    },
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

    let n1_hyper = n1.to_hypervolume_uv();
    let n2_hyper = n2.to_hypervolume_wx();

    let cross = n1_hyper.cross(&n2_hyper);

    // Check if the dot product
    let _hull = cross.get_convex_hull()?;

    Ok(vec![])
}

fn extend_if_not_present(result: &mut Vec<Point>, points: Vec<Point>) {
    for point in points.iter() {
        if !result.contains(point) {
            result.push(*point);
        } else {
            let index = result.iter().position(|p| *p == *point).unwrap();
            let existing_point = result.remove(index);
            let union = existing_point.union(*point);
            result.push(union);
        }
    }
}

struct RecursionCounter {
    count: usize,
}

fn refine(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
    resultant: &BernsteinHyperVolume<EFloat64>,
    split_u: bool,
    debug_data: BernsteinHyperVolume<Point>,
    counter: &mut RecursionCounter,
    scene_recorder: &mut PrimitiveSceneRecorder,
) -> AlgebraResult<Vec<Point>> {
    if counter.count == 0 {
        return Err(AlgebraError::new(
            "Recursion limit reached. Likely the solution is a curve".to_string(),
        ));
    }
    counter.count -= 1;

    let hull1 = s1.get_convex_hull();
    let hull2 = s2.get_convex_hull();

    match (hull1, hull2) {
        (Ok(hull1), Ok(hull2)) => {
            if counter.count > MAX_RECURSION_COUNT - 200 {
                let mut scene = PrimitiveScene::new();
                scene.add_surface_like(s1, Color10::Red, 10)?;
                scene.add_surface_like(s2, Color10::Blue, 10)?;

                scene.add_convex_hull(hull1.clone(), Color10::Green);
                scene.add_convex_hull(hull2.clone(), Color10::Green);

                if counter.count == MAX_RECURSION_COUNT - 20 {
                    scene.add_convex_hull(debug_data.get_convex_hull()?, Color10::Red);
                    scene.add_debug_text("debug_data".to_string());
                }

                for c in resultant.coefficients.iter() {
                    for c2 in c.iter() {
                        for c3 in c2.iter() {
                            for c4 in c3.iter() {
                                scene.add_point(
                                    Point::new(EFloat64::zero(), *c4, EFloat64::zero()),
                                    Color10::Green,
                                );
                            }
                        }
                    }
                }

                if hull1.intersects(&hull2) {
                    scene.add_debug_text("Intersects".to_string());
                } else {
                    scene.add_debug_text("Does not intersect".to_string());
                }

                scene_recorder.add_scene(scene);
            }

            if !hull1.intersects(&hull2) {
                return Ok(vec![]);
            }
        }
        _ => {
            let union = s1
                .eval(EFloat64::from(0.5), EFloat64::from(0.5))
                .union(s2.eval(EFloat64::from(0.5), EFloat64::from(0.5)));
            return Ok(vec![union]);
        }
    }

    if !resultant.could_be_zero() {
        return Ok(vec![]);
    }

    if split_u {
        let (s1a, s1b) = s1.split_u()?;
        let (s2a, s2b) = s2.split_u()?;

        let (a, b) = resultant.subdivide_u(EFloat64::from(0.5));
        let (raa, rab) = a.subdivide_w(EFloat64::from(0.5));
        let (rba, rbb) = b.subdivide_w(EFloat64::from(0.5));

        let (debug_a, debug_b) = debug_data.subdivide_u(EFloat64::from(0.5));
        let (debug_raa, debug_rab) = debug_a.subdivide_w(EFloat64::from(0.5));
        let (debug_rba, debug_rbb) = debug_b.subdivide_w(EFloat64::from(0.5));

        let mut result = vec![];

        extend_if_not_present(
            &mut result,
            refine(
                &s1a,
                &s2a,
                &raa,
                !split_u,
                debug_raa,
                counter,
                scene_recorder,
            )?,
        );
        extend_if_not_present(
            &mut result,
            refine(
                &s1a,
                &s2b,
                &rab,
                !split_u,
                debug_rab,
                counter,
                scene_recorder,
            )?,
        );
        extend_if_not_present(
            &mut result,
            refine(
                &s1b,
                &s2a,
                &rba,
                !split_u,
                debug_rba,
                counter,
                scene_recorder,
            )?,
        );
        extend_if_not_present(
            &mut result,
            refine(
                &s1b,
                &s2b,
                &rbb,
                !split_u,
                debug_rbb,
                counter,
                scene_recorder,
            )?,
        );

        return Ok(result);
    } else {
        let (s1a, s1b) = s1.split_v()?;
        let (s2a, s2b) = s2.split_v()?;

        let (a, b) = resultant.subdivide_v(EFloat64::from(0.5));
        let (raa, rab) = a.subdivide_x(EFloat64::from(0.5));
        let (rba, rbb) = b.subdivide_x(EFloat64::from(0.5));

        let (debug_a, debug_b) = debug_data.subdivide_v(EFloat64::from(0.5));
        let (debug_raa, debug_rab) = debug_a.subdivide_x(EFloat64::from(0.5));
        let (debug_rba, debug_rbb) = debug_b.subdivide_x(EFloat64::from(0.5));

        let mut result = vec![];

        extend_if_not_present(
            &mut result,
            refine(
                &s1a,
                &s2a,
                &raa,
                !split_u,
                debug_raa,
                counter,
                scene_recorder,
            )?,
        );
        extend_if_not_present(
            &mut result,
            refine(
                &s1a,
                &s2b,
                &rab,
                !split_u,
                debug_rab,
                counter,
                scene_recorder,
            )?,
        );
        extend_if_not_present(
            &mut result,
            refine(
                &s1b,
                &s2a,
                &rba,
                !split_u,
                debug_rba,
                counter,
                scene_recorder,
            )?,
        );
        extend_if_not_present(
            &mut result,
            refine(
                &s1b,
                &s2b,
                &rbb,
                !split_u,
                debug_rbb,
                counter,
                scene_recorder,
            )?,
        );

        return Ok(result);
    }
}

const MAX_RECURSION_COUNT: usize = 700;

pub fn get_critical_points_perpendicular_tangent_x(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
    scene_recorder: &mut PrimitiveSceneRecorder,
) -> AlgebraResult<Vec<Point>> {
    let hull_s1 = s1.get_convex_hull()?;
    let hull_s2 = s2.get_convex_hull()?;

    if !hull_s1.intersects(&hull_s2) {
        return Ok(vec![]);
    }

    // Now get normals
    let n1 = s1.normal();
    let n2 = s2.normal();

    let n1_hyper = n1.to_hypervolume_uv();
    let n2_hyper = n2.to_hypervolume_wx();

    let intersection_tangent = n1_hyper.cross(&n2_hyper);
    let product = intersection_tangent.dot(&BernsteinHyperVolume::unit_x());

    let result = refine(
        s1,
        s2,
        &product,
        true,
        intersection_tangent,
        &mut RecursionCounter {
            count: MAX_RECURSION_COUNT,
        },
        scene_recorder,
    )?;

    Ok(result)
}

pub fn get_critical_points_perpendicular_tangent_y(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
    scene_recorder: &mut PrimitiveSceneRecorder,
) -> AlgebraResult<Vec<Point>> {
    let hull_s1 = s1.get_convex_hull()?;
    let hull_s2 = s2.get_convex_hull()?;

    if !hull_s1.intersects(&hull_s2) {
        return Ok(vec![]);
    }

    // Now get normals
    let n1 = s1.normal();
    let n2 = s2.normal();

    let n1_hyper = n1.to_hypervolume_uv();
    let n2_hyper = n2.to_hypervolume_wx();

    let intersection_tangent = n1_hyper.cross(&n2_hyper);
    let product = intersection_tangent.dot(&BernsteinHyperVolume::unit_y());

    let result = refine(
        s1,
        s2,
        &product,
        true,
        intersection_tangent,
        &mut RecursionCounter {
            count: MAX_RECURSION_COUNT,
        },
        scene_recorder,
    )?;

    Ok(result)
}

pub fn get_critical_points_perpendicular_tangent_z(
    s1: &BernsteinSurface<Point>,
    s2: &BernsteinSurface<Point>,
    scene_recorder: &mut PrimitiveSceneRecorder,
) -> AlgebraResult<Vec<Point>> {
    let hull_s1 = s1.get_convex_hull()?;
    let hull_s2 = s2.get_convex_hull()?;

    if !hull_s1.intersects(&hull_s2) {
        return Ok(vec![]);
    }

    // Now get normals
    let n1 = s1.normal();
    let n2 = s2.normal();

    let n1_hyper = n1.to_hypervolume_uv();
    let n2_hyper = n2.to_hypervolume_wx();

    let intersection_tangent = n1_hyper.cross(&n2_hyper);
    let product = intersection_tangent.dot(&BernsteinHyperVolume::unit_z());

    let result = refine(
        s1,
        s2,
        &product,
        true,
        intersection_tangent,
        &mut RecursionCounter {
            count: MAX_RECURSION_COUNT,
        },
        scene_recorder,
    )?;

    Ok(result)
}

#[cfg(test)]
mod tests {
    use crate::{
        bernstein::bernstein_curve::BernsteinCurve,
        primitives::{color::Color10, primitive_scene::PrimitiveScene},
    };

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
    fn test_get_critical_points_perpendicular_tangent_x() -> AlgebraResult<()> {
        let s1 = mountain();
        let s2 = upside_down_mountain();

        let mut scene_recorder_x = PrimitiveSceneRecorder::new();
        let mut scene_recorder_y = PrimitiveSceneRecorder::new();

        let mut critical_points =
            get_critical_points_perpendicular_tangent_x(&s1, &s2, &mut scene_recorder_x)?;
        match get_critical_points_perpendicular_tangent_y(&s1, &s2, &mut scene_recorder_y) {
            Ok(points) => {
                critical_points.extend(points);
            }
            Err(_) => {
                println!("No critical points perpendicular to tangent y");
            }
        }
        // critical_points.extend(get_critical_points_perpendicular_tangent_z(&s1, &s2)?);

        println!("Critical points: {:?}", critical_points.len());
        scene_recorder_x.save_to_folder("test_outputs/critical_points_perpendicular_tangent_x")?;
        scene_recorder_y.save_to_folder("test_outputs/critical_points_perpendicular_tangent_y")?;

        let mut scene = PrimitiveScene::new();
        scene.add_surface_like(&s1, Color10::Red, 30)?;
        scene.add_surface_like(&s2, Color10::Blue, 30)?;

        for c in critical_points.iter() {
            scene.add_point(*c, Color10::Green);
        }

        scene.save_to_file("test_outputs/critical_points_perpendicular_tangent_x.html")?;

        Ok(())
    }

    #[test]
    fn test_get_distance_hypervolume() -> AlgebraResult<()> {
        let s1 = mountain();
        let s2 = upside_down_mountain();

        let distance = s1.to_hypervolume_uv() - s2.to_hypervolume_wx();
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

    #[test]
    fn subdivision_of_curve_example() -> AlgebraResult<()> {
        let s1 = BernsteinCurve::<Point>::new(vec![
            Point::from_f64(0.0, 0.0, 0.0),
            Point::from_f64(1.0, 3.0, 0.0),
            Point::from_f64(2.0, 1.0, 0.0),
        ]);

        let (mut a, mut b) = s1.subdivide(EFloat64::from(0.5));
        // shift a by 1 in z
        for c in a.coefficients.iter_mut() {
            c.z = c.z + EFloat64::one();
        }
        for c in b.coefficients.iter_mut() {
            c.z = c.z + EFloat64::one();
        }

        let mut scene = PrimitiveScene::new();
        scene.add_curve(&s1, Color10::Red)?;
        scene.add_curve(&a, Color10::Green)?;
        scene.add_curve(&b, Color10::Blue)?;

        for c in a.coefficients.iter() {
            scene.add_point(*c, Color10::Green);
        }
        for c in b.coefficients.iter() {
            scene.add_point(*c, Color10::Blue);
        }
        for c in s1.coefficients.iter() {
            scene.add_point(*c, Color10::Red);
        }

        scene.save_to_file("test_outputs/subdivision_of_curve_example.html")?;

        Ok(())
    }
}
