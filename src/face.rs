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
    use crate::{primitives::{color::{Color10, COLORS}, efloat::EFloat64, primitive_scene::PrimitiveScene}, surfaces::surface_like::SurfaceLike};

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
        let face = Face::try_new_from_surface(surface.clone())?;

        let mut scene = PrimitiveScene::new();
        scene.add_face(&face, Color10::Blue)?;
        scene.add_convex_hull(surface.get_convex_hull()?, Color10::Blue);

        scene.save_to_file("test_outputs/nurbs_surface.html")?;

        Ok(())
    }

    #[test]
    fn test_subdivide_surface() -> AlgebraResult<()> {
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

        let colors = COLORS;

        let mut faces = vec![surface];

        for level in 1..=3 {
            let mut new_faces = Vec::new();
            println!("Subdividing level {}", level);
            for surf in faces {
                let u_span = surf.u_span();
                let v_span = surf.v_span();
                let u_split = ((u_span.1 + u_span.0) / EFloat64::from(2.0))?;
                let v_split = ((v_span.1 + v_span.0) / EFloat64::from(2.0))?;
                let (surf1, surf2) = surf.subdivide_u(u_split)?;
                let (surf1a, surf1b) = surf1.subdivide_v(v_split)?;
                let (surf2a, surf2b) = surf2.subdivide_v(v_split)?;
                new_faces.push(surf1a);
                new_faces.push(surf1b);
                new_faces.push(surf2a);
                new_faces.push(surf2b);
            }

            println!("Rasterizing {} new faces", new_faces.len());
            faces = new_faces;
            let mut scene = PrimitiveScene::new();
            for (i, surf) in faces.iter().enumerate() {
                let face = Face::try_new_from_surface(surf.clone())?;
                scene.add_face(&face, colors[i % colors.len()])?;
                scene.add_convex_hull(surf.get_convex_hull()?, colors[i % colors.len()]);
            }
            let filename = format!("test_outputs/nurbs_surface_subdivide_level_{}.html", level);
            scene.save_to_file(&filename)?;
        }
        Ok(())
    }
}
