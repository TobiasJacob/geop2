use crate::{
    algebra_error::AlgebraResult,
    curves::curve_like::CurveLike,
    primitives::{efloat::EFloat64, line::Line},
};

pub fn rasterize_curve(curve: &dyn CurveLike) -> AlgebraResult<Vec<Line>> {
    let mut lines = Vec::new();
    let n = 30;
    let mut prev_point = curve.eval(EFloat64::from(0.0));
    for i in 0..n {
        let t = EFloat64::from(((i + 1) as f64 / n as f64).clamp(0.0, 1.0));
        let point = curve.eval(t);
        lines.push(Line::try_new(prev_point, point)?);
        prev_point = point;
    }
    Ok(lines)
}
