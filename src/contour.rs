use std::fmt::{self, Display};

use crate::{
    algebra_error::{AlgebraError, AlgebraResult, WithContext},
    curves::nurbs_curve::NurbsCurve,
};

pub struct Contour {
    pub bounds: Vec<NurbsCurve>,
}

impl Contour {
    pub fn try_new(curves: Vec<NurbsCurve>) -> AlgebraResult<Self> {
        let context = |err: AlgebraError| {
            let mut msg = "Failed to create contour".to_string();
            for curve in &curves {
                msg.push_str(&format!("\n{}", curve));
            }
            err.with_context(msg)
        };

        // Check if each end of the curve is the start of another curve.
        for i in 0..curves.len() {
            if curves[i].end_point() != curves[(i + 1) % curves.len()].start_point() {
                return Err(AlgebraError::new(format!(
                    "End point {} does not match start point {} at index {}",
                    curves[i].start_point(),
                    curves[(i + 1) % curves.len()].end_point(),
                    i
                )))
                .with_context(&context);
            }
        }
        Ok(Self { bounds: curves })
    }
}

impl Display for Contour {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Contour(")?;
        for bound in &self.bounds {
            write!(f, "{},", bound.start_point())?;
        }
        write!(f, ")")?;
        Ok(())
    }
}
