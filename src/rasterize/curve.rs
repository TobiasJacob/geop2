use crate::primitives::curve::CurveLike;

pub fn rasterize_curve(curve: &dyn CurveLike) -> Vec<Line> {
    let mut lines = Vec::new();
    let n = 30;
    for i in 0..=n {
        let t = EFloat64::from(i as f64 / n as f64);
        let point = curve.eval(t);
        lines.push(Line::new(point, point));
    }
    lines
}
