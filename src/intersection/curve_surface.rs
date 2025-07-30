use std::vec;

use crate::{
    algebra_error::AlgebraResult,
    curves::{curve_like::CurveLike, nurbs_curve::NurbsCurve},
    primitives::{
        color::Color10, efloat::EFloat64, point::Point, primitive_scene::PrimitiveScene,
        primitive_scene_recorder::PrimitiveSceneRecorder,
    },
    surfaces::{nurbs_surface::NurbsSurface, surface_like::SurfaceLike},
};

pub enum CurveCurveIntersection {
    Overlap(NurbsCurve),
    Points(Vec<Point>),
}

pub fn curve_surface_intersection_non_overlap(
    curve: &NurbsCurve,
    surface: &NurbsSurface,
    primitive_scene_recorder: &mut PrimitiveSceneRecorder,
) -> AlgebraResult<Vec<Point>> {
    // println!("--------------------------------");
    // println!("curve1: {}", curve1);
    // println!("curve2: {}", curve2);
    // let mut scene = PrimitiveScene::new();
    // scene.add_curve(curve1, Color10::Red)?;
    // scene.add_curve(curve2, Color10::Blue)?;
    // let cvx1 = curve1.control_polygon_hull()?;
    // let cvx2 = curve2.control_polygon_hull()?;
    // // if cvx1.is_err() {
    // //     primitive_scene_recorder.add_scene(scene);
    // //     return Ok(vec![curve1.coefficients[0]]);
    // // }
    // scene.add_convex_hull(cvx1, Color10::Red);
    // // if cvx2.is_err() {
    // //     primitive_scene_recorder.add_scene(scene);
    // //     return Ok(vec![curve2.coefficients[0]]);
    // // }
    // scene.add_convex_hull(cvx2, Color10::Blue);
    // primitive_scene_recorder.add_scene(scene);

    if !curve
        .control_polygon_hull()?
        .intersects(&surface.get_convex_hull()?)
    {
        return Ok(Vec::new());
    }

    let curve_span = curve.span();
    let surface_u_span = surface.u_span();
    let surface_v_span = surface.v_span();
    // println!("curve1_span: {:?}", curve1_span.1 - curve1_span.0);
    // println!("curve2_span: {:?}", curve2_span.1 - curve2_span.0);

    if curve_span.0 >= curve_span.1 - EFloat64::from(0.000001) {
        let p = curve.control_polygon_hull()?.to_point();
        println!("p: {:?}", p);
        return Ok(vec![p]);
    }

    if surface_u_span.0 >= surface_u_span.1 - EFloat64::from(0.000001)
        && surface_v_span.0 >= surface_v_span.1 - EFloat64::from(0.000001)
    {
        let p = surface.get_convex_hull()?.to_point();
        println!("p: {:?}", p);
        return Ok(vec![p]);
    }

    // If the split was not successful, we have reached the limit of the precision of the efloat
    if let Ok((curve1, curve2)) = curve.split() {
        if let Ok((surface1, surface2)) = surface.split_u() {
            if let Ok((surface1a, surface1b)) = surface1.split_v() {
                if let Ok((surface2a, surface2b)) = surface2.split_v() {
                    let mut scene = PrimitiveScene::new();
                    scene.add_curve(&curve1, Color10::Red).ok();
                    scene.add_curve(&curve2, Color10::Red).ok();
                    scene.add_surface(&surface1a, Color10::Blue, 10).ok();
                    scene.add_surface(&surface1b, Color10::Blue, 10).ok();
                    scene.add_surface(&surface2a, Color10::Blue, 10).ok();
                    scene.add_surface(&surface2b, Color10::Blue, 10).ok();

                    let curves = vec![("curve1", curve1), ("curve2", curve2)];
                    let surfaces = vec![
                        ("surface1a", surface1a),
                        ("surface1b", surface1b),
                        ("surface2a", surface2a),
                        ("surface2b", surface2b),
                    ];

                    let mut intersections = Vec::new();
                    for (name_curve, curve) in curves.iter() {
                        for (name_surface, surface) in surfaces.iter() {
                            if curve.control_polygon_hull()?.intersects(&surface.get_convex_hull()?) {
                                scene.add_convex_hull(curve.control_polygon_hull()?, Color10::Green);
                                scene.add_convex_hull(surface.get_convex_hull()?, Color10::Green);
                                scene.add_debug_text(format!(
                                    "{} intersects {}\n{:?}\n{:?}",
                                    name_curve, name_surface, curve.control_polygon_hull()?, surface.get_convex_hull()?
                                ));
                            }
                        }
                    }
                    primitive_scene_recorder.add_scene(scene);

                    for (_name_curve, curve) in curves.iter() {
                        for (_name_surface, surface) in surfaces.iter() {
                            intersections.extend(curve_surface_intersection_non_overlap(
                                &curve,
                                &surface,
                                primitive_scene_recorder,
                            )?);
                        }
                    }

                    // println!("curve1a: {}", curve1a);
                    // println!("curve1b: {}", curve1b);
                    // println!("curve2a: {}", curve2a);
                    // println!("curve2b: {}", curve2b);


                    // check if any of the intersections are the same point and if yes, unify them
                    let mut intersection_result = Vec::new();
                    for intersection in intersections {
                        let index = intersection_result
                            .iter()
                            .position(|p: &Point| *p == intersection);
                        match index {
                            Some(index) => {
                                println!(
                                    "union between: {} and {}",
                                    intersection_result[index], intersection
                                );
                                intersection_result[index] =
                                    intersection_result[index].union(intersection);
                            }
                            None => {
                                intersection_result.push(intersection);
                            }
                        }
                    }
                    return Ok(intersection_result);
                }
            }
        }
        let p = curve2.control_polygon_hull()?.to_point();
        println!("p: {:?}", p);
        return Ok(vec![p]);
    }
    let p = curve.control_polygon_hull()?.to_point();
    println!("p: {:?}", p);
    return Ok(vec![p]);
}

// test cases
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{intersection::test_surface_1, primitives::{color::Color10, efloat::EFloat64, primitive_scene::PrimitiveScene}};

    fn to_efloat_vec(values: Vec<f64>) -> Vec<EFloat64> {
        values.into_iter().map(EFloat64::from).collect()
    }

    #[test]
    fn test_curve_surface_intersection() -> AlgebraResult<()> {
        // Create control points
        let points1 = vec![
            Point::from_f64(0.0, 0.0, 0.0),
            Point::from_f64(1.0, 1.0, 0.0),
            Point::from_f64(2.0, 5.0, 0.0),
            Point::from_f64(3.0, 0.0, 0.0),
        ];

        // Create weights (all 1.0 for simplicity)
        let weights = to_efloat_vec(vec![1.0, 1.0, 1.0, 1.0]);

        // Create knot vector for degree 2 curve with 3 control points
        // For a degree 2 curve with 3 control points, we need 6 knots (n + p + 1 = 3 + 2 + 1 = 6)
        let knot_vector = to_efloat_vec(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]);

        // Create the curves
        let curve = NurbsCurve::try_new(points1, weights.clone(), knot_vector.clone(), 2)?;
        let surface = test_surface_1()?;

        let mut scene = PrimitiveScene::new();
        scene.add_curve(&curve, Color10::Red)?;
        scene.add_surface(&surface, Color10::Blue, 30)?;

        // Test intersection
        let mut primitive_scene_recorder = PrimitiveSceneRecorder::new();
        let points = curve_surface_intersection_non_overlap(&curve, &surface, &mut primitive_scene_recorder)?;

        for point in points {
            scene.add_point(point, Color10::Green);
            println!("point: {:?}", point);
        }

        primitive_scene_recorder.save_to_folder("test_outputs/curve_surface_intersection")?;
        scene.save_to_file("test_outputs/curve_surface_intersection.html")?;

        Ok(())
    }
}
