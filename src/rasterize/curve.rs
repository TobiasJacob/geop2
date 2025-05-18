use crate::{
    algebra_error::{AlgebraError, AlgebraResult, WithContext},
    curves::curve_like::CurveLike,
    primitives::{efloat::EFloat64, line::Line},
};

pub fn rasterize_curve(curve: &dyn CurveLike) -> AlgebraResult<Vec<Line>> {
    let context = |err: AlgebraError| err.with_context(format!("rasterizing curve: {}", curve));

    let span = curve.span().with_context(&context)?;
    let mut lines = Vec::new();
    let n = 30;
    let mut prev_point = curve.eval(span.0);
    for i in 0..n {
        let t = EFloat64::from(((i + 1) as f64 / n as f64).clamp(0.0, 1.0));
        let point = curve.eval(t * span.1 + (EFloat64::one() - t) * span.0);
        lines.push(Line::try_new(prev_point, point).with_context(&context)?);
        prev_point = point;
    }
    Ok(lines)
}
