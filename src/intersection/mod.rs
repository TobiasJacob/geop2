
pub mod curve_curve;
pub mod surface_surface;
pub mod curve_surface;

use crate::{algebra_error::AlgebraResult, primitives::{efloat::EFloat64, point::Point}, surfaces::nurbs_surface::NurbsSurface};
pub fn test_surface_1() -> AlgebraResult<NurbsSurface> {
    NurbsSurface::try_new(
        vec![
            vec![
                Point::from_f64(0.0, 2.0, 0.0),
                Point::from_f64(1.0, 1.0, 0.0),
                Point::from_f64(2.0, 0.0, 0.0),
                Point::from_f64(3.0, 1.0, 0.0),
            ],
            vec![
                Point::from_f64(0.0, 1.0, 1.0),
                Point::from_f64(1.0, 0.0, 1.0),
                Point::from_f64(2.0, 1.0, 1.0),
                Point::from_f64(3.0, 0.0, 1.0),
            ],
            vec![
                Point::from_f64(0.0, 0.0, 2.0),
                Point::from_f64(1.0, 1.0, 2.0),
                Point::from_f64(2.0, 0.0, 2.0),
                Point::from_f64(3.0, 1.0, 2.0),
            ],
            vec![
                Point::from_f64(0.0, 1.0, 3.0),
                Point::from_f64(1.0, 0.0, 3.0),
                Point::from_f64(2.0, 1.0, 3.0),
                Point::from_f64(3.0, 0.0, 3.0),
            ],
        ],
        vec![
            vec![
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
            ],
            vec![
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
            ],
            vec![
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
            ],
            vec![
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
                EFloat64::from(1.0),
            ],
        ],
        vec![
            EFloat64::from(0.0),
            EFloat64::from(0.0),
            EFloat64::from(0.0),
            EFloat64::from(0.5),
            EFloat64::from(1.0),
            EFloat64::from(1.0),
            EFloat64::from(1.0),
        ],
        vec![
            EFloat64::from(0.0),
            EFloat64::from(0.0),
            EFloat64::from(0.0),
            EFloat64::from(0.5),
            EFloat64::from(1.0),
            EFloat64::from(1.0),
            EFloat64::from(1.0),
        ],
        2,
        2,
    )
}
