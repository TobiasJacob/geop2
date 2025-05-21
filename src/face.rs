use std::fmt::{self, Display};

use crate::{
    algebra_error::{AlgebraError, AlgebraResult},
    contour::Contour,
    curves::nurbs_curve::NurbsCurve,
    primitives::point::Point,
    surfaces::nurbs_surface::NurbsSurface,
};

pub struct Face {
    pub surface: NurbsSurface,
    pub bounds: Vec<Contour>,
}

impl Face {
    pub fn try_new(surface: NurbsSurface, bounds: Vec<Contour>) -> AlgebraResult<Self> {
        // Check there is at least one contour
        if bounds.is_empty() {
            return Err(AlgebraError::new("No bounds found".to_string()));
        }

        Ok(Self { surface, bounds })
    }

    pub fn try_new_from_surface(surface: NurbsSurface) -> AlgebraResult<Self> {
        // Create a contour that goes from 0,0 to 1,0 to 1,1 to 0,1 to 0,0
        let bounds = Contour::try_new(vec![
            NurbsCurve::new_line(
                Point::from_f64(0.0, 0.0, 0.0),
                Point::from_f64(1.0, 0.0, 0.0),
            )?,
            NurbsCurve::new_line(
                Point::from_f64(1.0, 0.0, 0.0),
                Point::from_f64(1.0, 1.0, 0.0),
            )?,
            NurbsCurve::new_line(
                Point::from_f64(1.0, 1.0, 0.0),
                Point::from_f64(0.0, 1.0, 0.0),
            )?,
            NurbsCurve::new_line(
                Point::from_f64(0.0, 1.0, 0.0),
                Point::from_f64(0.0, 0.0, 0.0),
            )?,
        ])?;
        Self::try_new(surface, vec![bounds])
    }
}

impl Display for Face {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Face")?;
        write!(f, "Surface: {}", self.surface)?;
        write!(f, "Bounds: {}", self.bounds.len())?;
        for bound in &self.bounds {
            write!(f, "{}", bound)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::primitives::{color::Color10, efloat::EFloat64, primitive_scene::PrimitiveScene};

    use super::*;

    #[test]
    fn test_try_new_from_surface() -> AlgebraResult<()> {
        let surface = NurbsSurface::try_new(
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
        )?;
        let face = Face::try_new_from_surface(surface)?;

        let mut scene = PrimitiveScene::new();
        scene.add_face(&face, Color10::Blue)?;

        scene.save_to_file("test_outputs/nurbs_surface.html")?;

        Ok(())
    }
}
