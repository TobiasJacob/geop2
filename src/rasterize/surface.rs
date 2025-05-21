use crate::{
    algebra_error::{AlgebraError, AlgebraResult, WithContext},
    face::Face,
    primitives::{efloat::EFloat64, triangle::TriangleFace},
    surfaces::surface_like::SurfaceLike,
};

pub fn rasterize_surface(face: &Face) -> AlgebraResult<Vec<TriangleFace>> {
    let context =
        |err: AlgebraError| err.with_context(format!("rasterizing surface: {}", face.surface));

    let mut triangles = Vec::new();
    let n = 30;
    for i in 0..n {
        for j in 0..n {
            let point1 = face.surface.eval(
                EFloat64::from(i as f64 / n as f64),
                EFloat64::from(j as f64 / n as f64),
            );
            let point2 = face.surface.eval(
                EFloat64::from((i + 1) as f64 / n as f64),
                EFloat64::from(j as f64 / n as f64),
            );
            let point3 = face.surface.eval(
                EFloat64::from(i as f64 / n as f64),
                EFloat64::from((j + 1) as f64 / n as f64),
            );
            let point4 = face.surface.eval(
                EFloat64::from((i + 1) as f64 / n as f64),
                EFloat64::from((j + 1) as f64 / n as f64),
            );
            triangles.push(TriangleFace::try_new(point1, point2, point3).with_context(&context)?);
            triangles.push(TriangleFace::try_new(point3, point2, point4).with_context(&context)?);
        }
    }
    Ok(triangles)
}
