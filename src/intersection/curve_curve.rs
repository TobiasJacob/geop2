use crate::{
    algebra_error::AlgebraResult, curves::nurbs_curve::NurbsCurve, primitives::point::Point,
};

pub enum CurveCurveIntersection {
    Overlap(NurbsCurve),
    Point(Point),
}

pub fn curve_curve_intersection(
    curve1: NurbsCurve,
    curve2: NurbsCurve,
) -> AlgebraResult<CurveCurveIntersection> {
    // Check if the curves are the same

    todo!()
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
        scene.save_to_file("test_outputs/curve_curve_intersection.html")?;

        // Test intersection
        // let intersection = curve_curve_intersection(curve1, curve2)?;
        // assert!(matches!(intersection, CurveCurveIntersection::Overlap(_)));

        Ok(())
    }
}
