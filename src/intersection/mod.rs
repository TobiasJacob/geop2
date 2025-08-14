pub mod bernstein_surface;
pub mod curve_curve;
pub mod curve_surface;
pub mod surface_surface;

use crate::{
    algebra_error::AlgebraResult,
    primitives::{efloat::EFloat64, point::Point},
    surfaces::nurbs_surface::NurbsSurface,
};
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

#[cfg(test)]
mod tests {
    use crate::{
        bernstein::bernstein_curve::BernsteinCurve,
        primitives::{color::Color10, primitive_scene::PrimitiveScene},
    };

    use super::*;

    #[test]
    fn test_curve_curve_distance_surface() -> AlgebraResult<()> {
        let curve1 = BernsteinCurve::<Point>::new(vec![
            Point::from_f64(0.0, 0.0, 4.0),
            Point::from_f64(1.0, 1.0, 0.0),
            Point::from_f64(2.0, 5.0, 3.0),
            Point::from_f64(3.0, 1.0, 0.0),
        ]);
        let curve2 = BernsteinCurve::<Point>::new(vec![
            Point::from_f64(0.0, 1.0, 0.0),
            Point::from_f64(1.0, 1.0, 4.0),
            Point::from_f64(1.0, 1.0, 2.0),
            Point::from_f64(1.0, 1.0, 0.0),
        ]);

        let surface = curve1.to_surface_u() - curve2.to_surface_v();
        println!("surface: {}", surface);

        let mut scene = PrimitiveScene::new();
        scene.add_curve(&curve1, Color10::Blue)?;
        scene.add_curve(&curve2, Color10::Red)?;
        scene.add_surface_like(&surface, Color10::Green, 20)?;
        scene.save_to_file("test_outputs/curve_curve_distance_surface.html")?;

        Ok(())
    }

    #[test]
    fn test_curve_curve_distance_surface_with_one_self_intersecting_curve() -> AlgebraResult<()> {
        let curve1 = BernsteinCurve::<Point>::new(vec![
            Point::from_f64(0.0, 0.0, 0.0),
            Point::from_f64(1.0, 1.0, 0.0),
            Point::from_f64(2.0, 0.0, 0.0),
            Point::from_f64(2.0, -1.0, 0.0),
            Point::from_f64(0.0, 1.0, 0.0),
        ]);

        let surface = curve1.to_surface_u() - curve1.to_surface_v();
        println!("surface: {}", surface);

        let mut surface_z_strech = surface.clone();
        for (i, c) in surface_z_strech.coefficients.iter_mut().enumerate() {
            for c in c.iter_mut() {
                c.z = c.z + EFloat64::from(i as f64 * 0.3);
            }
        }
        println!("surface_z_strech: {}", surface_z_strech);

        let mut scene = PrimitiveScene::new();
        scene.add_curve(&curve1, Color10::Blue)?;
        scene.add_surface_like_wireframe(&surface, Color10::Blue, 20)?;
        scene.add_surface_like(&surface_z_strech, Color10::Blue, 20)?;
        let iso_u_0 = surface.iso_u_at(EFloat64::from(0.0));
        scene.add_curve(&iso_u_0, Color10::Red)?;
        let iso_v_0 = surface.iso_v_at(EFloat64::from(0.0));
        scene.add_curve(&iso_v_0, Color10::Green)?;
        let iso_u_1 = surface.iso_u_at(EFloat64::from(1.0));
        scene.add_curve(&iso_u_1, Color10::Red)?;
        let iso_v_1 = surface.iso_v_at(EFloat64::from(1.0));
        scene.add_curve(&iso_v_1, Color10::Green)?;
        scene.save_to_file(
            "test_outputs/curve_curve_distance_surface_with_one_self_intersecting_curve.html",
        )?;

        Ok(())
    }
}
