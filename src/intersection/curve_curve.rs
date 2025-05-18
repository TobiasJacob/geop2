use crate::{
    algebra_error::AlgebraResult,
    curves::{curve_like::CurveLike, nurbs_curve::NurbsCurve},
    primitives::{
        color::Color10, efloat::EFloat64, point::Point, primitive_scene::PrimitiveScene,
        primitive_scene_recorder::PrimitiveSceneRecorder,
    },
};

pub enum CurveCurveIntersection {
    Overlap(NurbsCurve),
    Points(Vec<Point>),
}

fn curve_curve_intersection_non_overlap(
    curve1: &NurbsCurve,
    curve2: &NurbsCurve,
    primitive_scene_recorder: &mut PrimitiveSceneRecorder,
) -> AlgebraResult<Vec<Point>> {
    // println!("--------------------------------");
    println!("curve1: {}", curve1);
    println!("curve2: {}", curve2);
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

    if !curve1
        .control_polygon_hull()?
        .intersects(&curve2.control_polygon_hull()?)
    {
        return Ok(Vec::new());
    }

    let curve1_span = curve1.span()?;
    let curve2_span = curve2.span()?;
    println!("curve1_span: {:?}", curve1_span.1 - curve1_span.0);
    println!("curve2_span: {:?}", curve2_span.1 - curve2_span.0);

    if curve1_span.0 >= curve1_span.1 - EFloat64::from(0.00001) {
        let p = curve1.control_polygon_hull()?.to_point();
        println!("p: {:?}", p);
        return Ok(vec![p]);
    }

    if curve2_span.0 >= curve2_span.1 - EFloat64::from(0.00001) {
        let p = curve2.control_polygon_hull()?.to_point();
        println!("p: {:?}", p);
        return Ok(vec![p]);
    }

    // If the split was not successful, we have reached the limit of the precision of the efloat
    if let Ok((curve1a, curve1b)) = curve1.split() {
        if let Ok((curve2a, curve2b)) = curve2.split() {
            let mut scene = PrimitiveScene::new();
            scene.add_curve(&curve1a, Color10::Red)?;
            scene.add_curve(&curve1b, Color10::Red)?;
            scene.add_curve(&curve2a, Color10::Blue)?;
            scene.add_curve(&curve2b, Color10::Blue)?;
            if curve1a
                .control_polygon_hull()?
                .intersects(&curve2a.control_polygon_hull()?)
            {
                scene.add_convex_hull(curve1a.control_polygon_hull()?, Color10::Green);
                scene.add_convex_hull(curve2a.control_polygon_hull()?, Color10::Green);
                scene.add_debug_text(format!(
                    "curve1a intersects curve2a\n{}\n{}",
                    curve1a.control_polygon_hull()?,
                    curve2a.control_polygon_hull()?
                ));
            }
            if curve1a
                .control_polygon_hull()?
                .intersects(&curve2b.control_polygon_hull()?)
            {
                scene.add_convex_hull(curve1a.control_polygon_hull()?, Color10::Olive);
                scene.add_convex_hull(curve2b.control_polygon_hull()?, Color10::Olive);
                scene.add_debug_text(format!(
                    "curve1a intersects curve2b\n{}\n{}",
                    curve1a.control_polygon_hull()?,
                    curve2b.control_polygon_hull()?
                ));
            }

            if curve1b
                .control_polygon_hull()?
                .intersects(&curve2a.control_polygon_hull()?)
            {
                scene.add_convex_hull(curve1b.control_polygon_hull()?, Color10::Cyan);
                scene.add_convex_hull(curve2a.control_polygon_hull()?, Color10::Cyan);
                scene.add_debug_text(format!(
                    "curve1b intersects curve2a\n{}\n{}",
                    curve1b.control_polygon_hull()?,
                    curve2a.control_polygon_hull()?
                ));
            }
            if curve1b
                .control_polygon_hull()?
                .intersects(&curve2b.control_polygon_hull()?)
            {
                scene.add_convex_hull(curve1b.control_polygon_hull()?, Color10::Purple);
                scene.add_convex_hull(curve2b.control_polygon_hull()?, Color10::Purple);
                scene.add_debug_text(format!(
                    "curve1b intersects curve2b\n{}\n{}",
                    curve1b.control_polygon_hull()?,
                    curve2b.control_polygon_hull()?
                ));
            }
            primitive_scene_recorder.add_scene(scene);
            // println!("curve1a: {}", curve1a);
            // println!("curve1b: {}", curve1b);
            // println!("curve2a: {}", curve2a);
            // println!("curve2b: {}", curve2b);

            let mut intersections = Vec::new();
            intersections.extend(curve_curve_intersection_non_overlap(
                &curve1a,
                &curve2a,
                primitive_scene_recorder,
            )?);
            intersections.extend(curve_curve_intersection_non_overlap(
                &curve1a,
                &curve2b,
                primitive_scene_recorder,
            )?);
            intersections.extend(curve_curve_intersection_non_overlap(
                &curve1b,
                &curve2a,
                primitive_scene_recorder,
            )?);
            intersections.extend(curve_curve_intersection_non_overlap(
                &curve1b,
                &curve2b,
                primitive_scene_recorder,
            )?);

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
                        intersection_result[index] = intersection_result[index].union(intersection);
                    }
                    None => {
                        intersection_result.push(intersection);
                    }
                }
            }
            return Ok(intersection_result);
        }
        let p = curve2.control_polygon_hull()?.to_point();
        println!("p: {:?}", p);
        return Ok(vec![p]);
    }
    let p = curve1.control_polygon_hull()?.to_point();
    println!("p: {:?}", p);
    return Ok(vec![p]);
}

pub fn curve_curve_intersection(
    curve1: NurbsCurve,
    curve2: NurbsCurve,
) -> AlgebraResult<CurveCurveIntersection> {
    let mut primitive_scene_recorder = PrimitiveSceneRecorder::new();
    let intersections =
        curve_curve_intersection_non_overlap(&curve1, &curve2, &mut primitive_scene_recorder);
    primitive_scene_recorder.save_to_folder("test_outputs/curve_curve_intersection")?;
    let intersections = intersections?;
    Ok(CurveCurveIntersection::Points(intersections))
}

// test cases
#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{color::Color10, efloat::EFloat64, primitive_scene::PrimitiveScene};

    fn to_efloat_vec(values: Vec<f64>) -> Vec<EFloat64> {
        values.into_iter().map(EFloat64::from).collect()
    }

    #[test]
    fn test_curve_curve_intersection() -> AlgebraResult<()> {
        // Create control points
        let points1 = vec![
            Point::from_f64(0.0, 0.0, 0.0),
            Point::from_f64(1.0, 1.0, 0.0),
            Point::from_f64(2.0, 5.0, 0.0),
            Point::from_f64(3.0, 0.0, 0.0),
        ];

        let points2 = vec![
            Point::from_f64(0.0, 1.5, 0.0),
            Point::from_f64(3.0, 4.0, 0.0),
            Point::from_f64(2.0, 2.0, 0.0),
            Point::from_f64(3.0, 2.0, 0.0),
        ];

        // Create weights (all 1.0 for simplicity)
        let weights = to_efloat_vec(vec![1.0, 1.0, 1.0, 1.0]);

        // Create knot vector for degree 2 curve with 3 control points
        // For a degree 2 curve with 3 control points, we need 6 knots (n + p + 1 = 3 + 2 + 1 = 6)
        let knot_vector = to_efloat_vec(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]);

        // Create the curves
        let curve1 = NurbsCurve::try_new(points1, weights.clone(), knot_vector.clone(), 2)?;
        let curve2 = NurbsCurve::try_new(points2, weights, knot_vector, 2)?;

        let mut scene = PrimitiveScene::new();
        scene.add_curve(&curve1, Color10::Red)?;
        scene.add_curve(&curve2, Color10::Blue)?;

        // Test intersection
        let intersection = curve_curve_intersection(curve1, curve2)?;

        match intersection {
            CurveCurveIntersection::Overlap(_curve) => {}
            CurveCurveIntersection::Points(points) => {
                for point in points {
                    scene.add_point(point, Color10::Green);
                    println!("point: {:?}", point);
                }
            }
        }

        scene.save_to_file("test_outputs/curve_curve_intersection.html")?;

        Ok(())
    }
}
