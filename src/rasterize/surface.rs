use crate::{
    algebra_error::{AlgebraError, AlgebraResult, WithContext},
    face::Face,
    primitives::{efloat::EFloat64, triangle::TriangleFace},
    surfaces::surface_like::SurfaceLike,
};

pub fn rasterize_surface_like(
    surface: &dyn SurfaceLike,
    n: usize,
) -> AlgebraResult<Vec<TriangleFace>> {
    let context = |err: AlgebraError| err.with_context(format!("rasterizing surface: {}", surface));
    let mut triangles = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let point1 = surface.eval(
                EFloat64::from(i as f64 / n as f64),
                EFloat64::from(j as f64 / n as f64),
            );
            let point2 = surface.eval(
                EFloat64::from((i + 1) as f64 / n as f64),
                EFloat64::from(j as f64 / n as f64),
            );
            let point3 = surface.eval(
                EFloat64::from(i as f64 / n as f64),
                EFloat64::from((j + 1) as f64 / n as f64),
            );
            let point4 = surface.eval(
                EFloat64::from((i + 1) as f64 / n as f64),
                EFloat64::from((j + 1) as f64 / n as f64),
            );
            let t1 = TriangleFace::try_new(point1, point2, point3).with_context(&context);
            if let Ok(t1) = t1 {
                triangles.push(t1);
            }
            let t2 = TriangleFace::try_new(point3, point2, point4).with_context(&context);
            if let Ok(t2) = t2 {
                triangles.push(t2);
            }
        }
    }
    Ok(triangles)
}

pub fn rasterize_face(face: &Face, n: usize) -> AlgebraResult<Vec<TriangleFace>> {
    let context =
        |err: AlgebraError| err.with_context(format!("rasterizing surface: {}", face.surface));

    let mut triangles = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let span_u = face.surface.u_span();
            let diff_u = span_u.1 - span_u.0;
            let span_v = face.surface.v_span();
            let diff_v = span_v.1 - span_v.0;
            let point1 = face.surface.eval(
                EFloat64::from(i as f64 / n as f64) * diff_u + span_u.0,
                EFloat64::from(j as f64 / n as f64) * diff_v + span_v.0,
            );
            let point2 = face.surface.eval(
                EFloat64::from((i + 1) as f64 / n as f64) * diff_u + span_u.0,
                EFloat64::from(j as f64 / n as f64) * diff_v + span_v.0,
            );
            let point3 = face.surface.eval(
                EFloat64::from(i as f64 / n as f64) * diff_u + span_u.0,
                EFloat64::from((j + 1) as f64 / n as f64) * diff_v + span_v.0,
            );
            let point4 = face.surface.eval(
                EFloat64::from((i + 1) as f64 / n as f64) * diff_u + span_u.0,
                EFloat64::from((j + 1) as f64 / n as f64) * diff_v + span_v.0,
            );
            triangles.push(TriangleFace::try_new(point1, point2, point3).with_context(&context)?);
            triangles.push(TriangleFace::try_new(point3, point2, point4).with_context(&context)?);
        }
    }
    Ok(triangles)
}
