use std::fmt::Display;

use crate::{
    algebra_error::{AlgebraError, AlgebraResult},
    primitives::{efloat::EFloat64, point::Point},
};

/// A B‑spline surface defined as a tensor product of B‑spline basis functions.
///
/// The control net is stored as a Vec<Vec<Point>> (rows × columns). There are two knot vectors,
/// one for the u–direction and one for the v–direction, with corresponding degrees.
#[derive(Debug, Clone)]
pub struct BSplineSurface {
    pub coefficients: Vec<Vec<Point>>, // 2D grid of control points.
    pub knot_vector_u: Vec<EFloat64>,  // Knot vector in the u-direction.
    pub knot_vector_v: Vec<EFloat64>,  // Knot vector in the v-direction.
    pub degree_u: usize,               // Degree in the u-direction.
    pub degree_v: usize,               // Degree in the v-direction.
}

impl BSplineSurface {
    /// Constructs a new BSplineSurface.
    ///
    /// # Parameters
    /// - `coefficients`: A 2D grid (rows × columns) of control points.
    /// - `knot_vector_u`: Knot vector for the u–direction.
    /// - `knot_vector_v`: Knot vector for the v–direction.
    /// - `degree_u`: Degree for the u–direction.
    /// - `degree_v`: Degree for the v–direction.
    ///
    /// # Errors
    /// Returns an error if the control net dimensions are insufficient for the given degrees or if
    /// the knot vectors are not of the correct length or not non-decreasing.
    pub fn try_new(
        coefficients: Vec<Vec<Point>>,
        knot_vector_u: Vec<EFloat64>,
        knot_vector_v: Vec<EFloat64>,
        degree_u: usize,
        degree_v: usize,
    ) -> AlgebraResult<Self> {
        if coefficients.is_empty() {
            return Err(AlgebraError::new(
                "BSplineSurface invalid input: coefficients cannot be empty".to_string(),
            ));
        }
        let num_rows = coefficients.len();
        let num_cols = coefficients[0].len();
        // Ensure all rows have the same number of columns.
        for row in &coefficients {
            if row.len() != num_cols {
                return Err(AlgebraError::new(
                    "BSplineSurface invalid input: inconsistent number of columns in control points"
                        .to_string(),
                ));
            }
        }

        if num_rows <= degree_u {
            return Err(AlgebraError::new(format!(
                "BSplineSurface invalid input: number of rows ({}) <= degree_u ({})",
                num_rows, degree_u
            )));
        }
        if num_cols <= degree_v {
            return Err(AlgebraError::new(format!(
                "BSplineSurface invalid input: number of columns ({}) <= degree_v ({})",
                num_cols, degree_v
            )));
        }

        if knot_vector_u.len() != num_rows + 1 + degree_u {
            return Err(AlgebraError::new(format!(
                "BSplineSurface invalid input: knot_vector_u.len() ({}) != number of rows ({}) + 1 + degree_u ({})",
                knot_vector_u.len(),
                num_rows,
                degree_u
            )));
        }
        if knot_vector_v.len() != num_cols + 1 + degree_v {
            return Err(AlgebraError::new(format!(
                "BSplineSurface invalid input: knot_vector_v.len() ({}) != number of columns ({}) + 1 + degree_v ({})",
                knot_vector_v.len(),
                num_cols,
                degree_v
            )));
        }

        // Check that knot vectors are non-decreasing.
        for i in 1..knot_vector_u.len() {
            if knot_vector_u[i - 1] > knot_vector_u[i] {
                return Err("Knot vector u must be non-decreasing".into());
            }
        }
        for i in 1..knot_vector_v.len() {
            if knot_vector_v[i - 1] > knot_vector_v[i] {
                return Err("Knot vector v must be non-decreasing".into());
            }
        }

        Ok(Self {
            coefficients,
            knot_vector_u,
            knot_vector_v,
            degree_u,
            degree_v,
        })
    }

    /// Constructs a BSplineSurface corresponding to a tensor–product basis function.
    /// All control points are set to zero except for the one at (index_u, index_v) which is set
    /// to `unit_vector`.
    pub fn try_new_from_basis(
        index_u: usize,
        index_v: usize,
        degree_u: usize,
        degree_v: usize,
        knot_vector_u: Vec<EFloat64>,
        knot_vector_v: Vec<EFloat64>,
        unit_vector: Point,
    ) -> AlgebraResult<BSplineSurface> {
        let num_rows = knot_vector_u.len() - 1 - degree_u;
        let num_cols = knot_vector_v.len() - 1 - degree_v;
        if index_u >= num_rows {
            return Err(AlgebraError::new(format!(
                "BSplineSurface invalid input: index_u {} is out of bounds (num_rows = {})",
                index_u, num_rows
            )));
        }
        if index_v >= num_cols {
            return Err(AlgebraError::new(format!(
                "BSplineSurface invalid input: index_v {} is out of bounds (num_cols = {})",
                index_v, num_cols
            )));
        }
        let mut coefficients = vec![vec![Point::zero(); num_cols]; num_rows];
        coefficients[index_u][index_v] = unit_vector;
        Self::try_new(
            coefficients,
            knot_vector_u,
            knot_vector_v,
            degree_u,
            degree_v,
        )
    }

    /// A generic knot–span finder.
    /// Returns an index k such that knot_vector[k] ≤ t < knot_vector[k+1].
    fn find_span_generic(knot_vector: &[EFloat64], t: EFloat64) -> Option<usize> {
        if t < knot_vector[0] {
            return None;
        }
        if t >= knot_vector[knot_vector.len() - 1] {
            return None;
        }
        let mut span = 0;
        while !(knot_vector[span] <= t && t < knot_vector[span + 1]) {
            span += 1;
        }
        Some(span)
    }

    pub fn find_span_u(&self, t: EFloat64) -> Option<usize> {
        BSplineSurface::find_span_generic(&self.knot_vector_u, t)
    }

    pub fn find_span_v(&self, t: EFloat64) -> Option<usize> {
        BSplineSurface::find_span_generic(&self.knot_vector_v, t)
    }

    /// Evaluates the B‑spline surface at the parameter pair (u, v).
    ///
    /// This is done using a tensor product de Boor algorithm:
    /// 1. For each contributing row (in the u–direction), the surface is evaluated in the v–direction.
    /// 2. The resulting intermediate points are then used in a de Boor evaluation in the u–direction.
    pub fn eval(&self, u: EFloat64, v: EFloat64) -> Point {
        let k_u = match self.find_span_u(u.clone()) {
            Some(idx) => idx,
            None => return Point::zero(),
        };
        let k_v = match self.find_span_v(v.clone()) {
            Some(idx) => idx,
            None => return Point::zero(),
        };
        let p_u = self.degree_u;
        let p_v = self.degree_v;

        let mut d: Vec<Vec<Point>> = Vec::with_capacity(p_u + 1);
        for i in 0..=p_u {
            if k_u + i < p_u || k_u + i - p_u >= self.coefficients.len() {
                d.push(vec![Point::zero(); self.coefficients[0].len()]);
            } else {
                d.push(self.coefficients[k_u + i - p_u].clone());
                for j in 0..=p_v {
                    if k_v + j < p_v || k_v + j - p_v >= self.coefficients[0].len() {
                        d[i][j] = Point::zero();
                    } else {
                        d[i][j] = self.coefficients[k_u + i - p_u][k_v + j - p_v].clone();
                    }
                }
            }
        }

        // Apply de Boor's algorithm.
        for r_u in 1..=p_u {
            for j_u in (r_u..=p_u).rev() {
                let alpha_u =
                    match k_u + j_u < p_u || j_u + 1 + k_u - r_u >= self.knot_vector_u.len() {
                        true => EFloat64::zero(),
                        false => {
                            let left_knot = self.knot_vector_u[j_u + k_u - p_u].clone();
                            let right_knot = self.knot_vector_u[j_u + 1 + k_u - r_u].clone();

                            // In case we divide by zero, the alpha value does not matter, so we choose 0.
                            ((u - left_knot) / (right_knot - left_knot)).unwrap_or(EFloat64::zero())
                        }
                    };
                for j_v in 0..=p_v {
                    d[j_u][j_v] = d[j_u - 1][j_v].clone() * (EFloat64::one() - alpha_u.clone())
                        + d[j_u][j_v].clone() * alpha_u.clone();
                }
            }
        }

        for r_v in 1..=p_v {
            for j_v in (r_v..=p_v).rev() {
                let alpha_v =
                    match k_v + j_v < p_v || j_v + 1 + k_v - r_v >= self.knot_vector_v.len() {
                        true => EFloat64::zero(),
                        false => {
                            let left_knot = self.knot_vector_v[j_v + k_v - p_v].clone();
                            let right_knot = self.knot_vector_v[j_v + 1 + k_v - r_v].clone();

                            // In case we divide by zero, the alpha value does not matter, so we choose 0.
                            ((v - left_knot) / (right_knot - left_knot)).unwrap_or(EFloat64::zero())
                        }
                    };
                d[p_u][j_v] = d[p_u][j_v - 1].clone() * (EFloat64::one() - alpha_v.clone())
                    + d[p_u][j_v].clone() * alpha_v.clone();
            }
        }

        d[p_u][p_v].clone()
    }

    /// Inserts the knot value `t` once into the B‑spline surface in the u–direction using the standard
    /// knot insertion algorithm. This updates the u knot vector and control net.
    pub fn insert_knot_u(&self, t: EFloat64) -> AlgebraResult<BSplineSurface> {
        let k = match self.find_span_u(t.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Parameter t is out of the valid domain for knot insertion in u".to_string(),
                ));
            }
        };
        println!("insert_knot_u: k = {}", k);
        let p = self.degree_u;
        let m = self.coefficients.len() - 1; // last row index
        let num_cols = self.coefficients[0].len();

        // Build new u knot vector.
        let mut new_knot_vector_u: Vec<EFloat64> = Vec::with_capacity(self.knot_vector_u.len() + 1);
        for i in 0..=k {
            new_knot_vector_u.push(self.knot_vector_u[i].clone());
        }
        new_knot_vector_u.push(t.clone());
        for i in (k + 1)..self.knot_vector_u.len() {
            new_knot_vector_u.push(self.knot_vector_u[i].clone());
        }

        // Build new control net (u–direction only changes: columns remain unchanged).
        let mut new_coeffs: Vec<Vec<Point>> = Vec::with_capacity(self.coefficients.len() + 1);
        // Copy rows that are not affected.
        for i in 0..(k - p + 1) {
            new_coeffs.push(self.coefficients[i].clone());
        }
        // Recompute affected rows.
        for i in (k - p + 1)..=k {
            let mut new_row: Vec<Point> = Vec::with_capacity(num_cols);
            let alpha = ((t.clone() - self.knot_vector_u[i].clone())
                / (self.knot_vector_u[i + p].clone() - self.knot_vector_u[i].clone()))
            .unwrap_or(EFloat64::zero());
            for j in 0..num_cols {
                let pt = self.coefficients[i - 1][j].clone() * (EFloat64::one() - alpha.clone())
                    + self.coefficients[i][j].clone() * alpha.clone();
                new_row.push(pt);
            }
            new_coeffs.push(new_row);
        }
        // Copy the remaining rows.
        for i in k..=m {
            new_coeffs.push(self.coefficients[i].clone());
        }

        BSplineSurface::try_new(
            new_coeffs,
            new_knot_vector_u,
            self.knot_vector_v.clone(),
            p,
            self.degree_v,
        )
    }

    /// Subdivides the B‑spline surface along the u–direction at parameter `t` into two new
    /// BSplineSurface segments. The method inserts `t` repeatedly until its multiplicity equals
    /// degree_u+1 (i.e. a break point) then splits the control net and u knot vector.
    ///
    /// The v–direction remains unchanged.
    pub fn subdivide_u(&self, t: EFloat64) -> AlgebraResult<(BSplineSurface, BSplineSurface)> {
        let p = self.degree_u;

        if t < self.knot_vector_u[p] || t > self.knot_vector_u[self.knot_vector_u.len() - p - 1] {
            return Err(AlgebraError::new(
                "Parameter t is out of the valid domain for u subdivision".to_string(),
            ));
        }

        // Determine current multiplicity of t in u–direction.
        let current_multiplicity = self.knot_vector_u.iter().filter(|&knot| *knot == t).count();
        // t must appear with multiplicity p+1.
        let r = p - current_multiplicity;
        let mut surface = self.clone();
        for _ in 0..r {
            println!("subdivide_u: surface = {}", surface);
            surface = surface.insert_knot_u(t.clone())?;
        }

        let t_index = match surface.find_span_u(t.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Failed to locate knot t with full multiplicity after insertion in u"
                        .to_string(),
                ));
            }
        };

        // Build left surface: rows 0 to (t_index - p) and the corresponding u–knot vector.
        let left_coeffs = surface.coefficients[..=(t_index - p)].to_vec();
        let left_knots_u = surface.knot_vector_u[..=t_index + 1].to_vec();
        let left_surface = BSplineSurface::try_new(
            left_coeffs,
            left_knots_u,
            surface.knot_vector_v.clone(),
            p,
            surface.degree_v,
        )?;

        // Build right surface: rows from (t_index - p) to end.
        let right_coeffs = surface.coefficients[(t_index - p)..].to_vec();
        let right_knots_u = surface.knot_vector_u[t_index - p..].to_vec();
        let right_surface = BSplineSurface::try_new(
            right_coeffs,
            right_knots_u,
            surface.knot_vector_v.clone(),
            p,
            surface.degree_v,
        )?;

        Ok((left_surface, right_surface))
    }
}

impl Display for BSplineSurface {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "BSplineSurface(")?;
        writeln!(f, "  Control Points:")?;
        for row in &self.coefficients {
            write!(f, "    [")?;
            for pt in row {
                write!(f, "{}, ", pt)?;
            }
            writeln!(f, "]")?;
        }
        writeln!(
            f,
            "  Degrees: (u: {}, v: {}),",
            self.degree_u, self.degree_v
        )?;
        write!(f, "  Knot Vector u: [")?;
        for knot in &self.knot_vector_u {
            write!(f, "{}, ", knot)?;
        }
        write!(f, "],\n  Knot Vector v: [")?;
        for knot in &self.knot_vector_v {
            write!(f, "{}, ", knot)?;
        }
        write!(f, "])")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_efloat_vec(values: Vec<f64>) -> Vec<EFloat64> {
        values.into_iter().map(EFloat64::from).collect()
    }

    #[test]
    fn test_eval() -> AlgebraResult<()> {
        // Create a simple surface (3 rows × 3 columns) with control points along unit z scaled by a factor.
        let coefficients = vec![
            vec![
                Point::unit_z() * EFloat64::from(1.0),
                Point::unit_z() * EFloat64::from(2.0),
                Point::unit_z() * EFloat64::from(3.0),
            ],
            vec![
                Point::unit_z() * EFloat64::from(4.0),
                Point::unit_z() * EFloat64::from(5.0),
                Point::unit_z() * EFloat64::from(6.0),
            ],
            vec![
                Point::unit_z() * EFloat64::from(7.0),
                Point::unit_z() * EFloat64::from(8.0),
                Point::unit_z() * EFloat64::from(9.0),
            ],
        ];

        // Use clamped knot vectors.
        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface = BSplineSurface::try_new(coefficients, knot_vector_u, knot_vector_v, 2, 2)?;

        println!("Surface: {}", surface);

        let u = EFloat64::from(0.5);
        let v = EFloat64::from(0.5);
        let result = surface.eval(u, v);
        println!("Surface eval at (0.5, 0.5): {}", result);
        assert!(result != Point::zero());

        // test at 0 and 1
        let u = EFloat64::from(0.0);
        let v = EFloat64::from(0.0);
        let result = surface.eval(u, v);
        println!("Surface eval at (0.0, 0.0): {}", result);
        assert_eq!(result, Point::unit_z() * EFloat64::from(1.0));

        let u = EFloat64::from(0.9999);
        let v = EFloat64::from(0.0);
        let result = surface.eval(u, v);
        println!("Surface eval at (1.0, 0.0): {}", result);
        assert_eq!(result, Point::unit_z() * EFloat64::new(7.0, 6.999));

        let u = EFloat64::from(0.0);
        let v = EFloat64::from(0.9999);
        let result = surface.eval(u, v);
        println!("Surface eval at (0.0, 1.0): {}", result);
        assert_eq!(result, Point::unit_z() * EFloat64::new(3.0, 2.999));

        let u = EFloat64::from(0.9999);
        let v = EFloat64::from(0.9999);
        let result = surface.eval(u, v);
        println!("Surface eval at (1.0, 1.0): {}", result);
        assert_eq!(result, Point::unit_z() * EFloat64::new(9.0, 8.999));

        Ok(())
    }

    #[test]
    fn test_insert_knot_u() -> AlgebraResult<()> {
        // Create a surface similar to test_eval.
        let coefficients = vec![
            vec![
                Point::unit_z() * EFloat64::from(1.0),
                Point::unit_z() * EFloat64::from(2.0),
                Point::unit_z() * EFloat64::from(3.0),
            ],
            vec![
                Point::unit_z() * EFloat64::from(4.0),
                Point::unit_z() * EFloat64::from(5.0),
                Point::unit_z() * EFloat64::from(6.0),
            ],
            vec![
                Point::unit_z() * EFloat64::from(7.0),
                Point::unit_z() * EFloat64::from(8.0),
                Point::unit_z() * EFloat64::from(9.0),
            ],
        ];
        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface = BSplineSurface::try_new(coefficients, knot_vector_u, knot_vector_v, 2, 2)?;
        println!("Original surface: {}", surface);
        let t = EFloat64::from(1.0);
        let surface2 = surface.insert_knot_u(t.clone())?;
        println!("Surface after inserting knot u at {}: {}", t, surface2);

        // Check that evaluation remains the same after knot insertion.
        for i in 0..=100 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            for j in 0..=100 {
                let v_val = EFloat64::from(j as f64 / 100.0);
                let orig = surface.eval(u_val.clone(), v_val.clone());
                let new = surface2.eval(u_val.clone(), v_val.clone());
                assert_eq!(orig, new);
            }
        }
        Ok(())
    }

    #[test]
    fn test_subdivide_u() -> AlgebraResult<()> {
        // Create a simple surface.
        let coefficients = vec![
            vec![
                Point::unit_z() * EFloat64::from(1.0),
                Point::unit_z() * EFloat64::from(2.0),
                Point::unit_z() * EFloat64::from(3.0),
            ],
            vec![
                Point::unit_z() * EFloat64::from(4.0),
                Point::unit_z() * EFloat64::from(5.0),
                Point::unit_z() * EFloat64::from(6.0),
            ],
            vec![
                Point::unit_z() * EFloat64::from(7.0),
                Point::unit_z() * EFloat64::from(8.0),
                Point::unit_z() * EFloat64::from(9.0),
            ],
        ];
        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface = BSplineSurface::try_new(coefficients, knot_vector_u, knot_vector_v, 2, 2)?;
        println!("Original surface: {}", surface);
        let t = EFloat64::from(1.0);
        let (left, right) = surface.subdivide_u(t.clone())?;
        println!("Left surface: {}", left);
        println!("Right surface: {}", right);

        // Check that evaluation on the left and right segments matches the original surface.
        for i in 0..=50 {
            let u_val = EFloat64::from(i as f64 / 50.0 * 1.0);
            for j in 0..=50 {
                let v_val = EFloat64::from(j as f64 / 50.0);
                let orig = surface.eval(u_val.clone(), v_val.clone());
                let left_val = left.eval(u_val.clone(), v_val.clone());
                assert_eq!(orig, left_val, "Mismatch at (u,v)=({}, {})", u_val, v_val);
            }
        }
        for i in 0..=50 {
            let u_val = EFloat64::from(i as f64 / 50.0 * 1.0 + 1.0);
            for j in 0..=50 {
                let v_val = EFloat64::from(j as f64 / 50.0);
                let orig = surface.eval(u_val.clone(), v_val.clone());
                let right_val = right.eval(u_val.clone(), v_val.clone());
                assert_eq!(orig, right_val, "Mismatch at (u,v)=({}, {})", u_val, v_val);
            }
        }

        Ok(())
    }
}
