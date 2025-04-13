use std::fmt::Display;

use crate::{
    algebra_error::{AlgebraError, AlgebraResult},
    primitives::{efloat::EFloat64, point::Point},
};

/// A NURBS (Non-Uniform Rational B-Spline) surface defined as a tensor product of NURBS basis functions.
///
/// - `coefficients` is a 2D grid (rows × columns) of control points.
/// - `weights` is a 2D grid of corresponding weights.
/// - `knot_vector_u` and `knot_vector_v` are the knot vectors for the u– and v–directions respectively.
/// - `degree_u` and `degree_v` are the polynomial degrees in the corresponding directions.
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    pub coefficients: Vec<Vec<Point>>, // 2D grid of control points.
    pub weights: Vec<Vec<EFloat64>>,   // 2D grid of weights corresponding to each control point.
    pub knot_vector_u: Vec<EFloat64>,  // Knot vector in the u–direction.
    pub knot_vector_v: Vec<EFloat64>,  // Knot vector in the v–direction.
    pub degree_u: usize,               // Degree in the u–direction.
    pub degree_v: usize,               // Degree in the v–direction.
}

impl NurbsSurface {
    /// Constructs a new NurbsSurface.
    ///
    /// Checks that:
    /// - The control net (coefficients) is nonempty and all rows have the same number of columns.
    /// - The weight matrix has the same dimensions as the control net.
    /// - There are enough rows/columns given the degrees.
    /// - The u– and v–knot vector lengths equal (rows/cols + 1 + degree) respectively and are non-decreasing.
    pub fn try_new(
        coefficients: Vec<Vec<Point>>,
        weights: Vec<Vec<EFloat64>>,
        knot_vector_u: Vec<EFloat64>,
        knot_vector_v: Vec<EFloat64>,
        degree_u: usize,
        degree_v: usize,
    ) -> AlgebraResult<Self> {
        if coefficients.is_empty() {
            return Err(AlgebraError::new(
                "NurbsSurface invalid input: coefficients cannot be empty".to_string(),
            ));
        }
        let num_rows = coefficients.len();
        let num_cols = coefficients[0].len();
        // Check that every row in coefficients and weights has the same number of columns.
        for (i, row) in coefficients.iter().enumerate() {
            if row.len() != num_cols {
                return Err(AlgebraError::new(
                    "NurbsSurface invalid input: inconsistent number of columns in coefficients"
                        .to_string(),
                ));
            }
            if weights.get(i).map_or(true, |w_row| w_row.len() != num_cols) {
                return Err(AlgebraError::new(
                    "NurbsSurface invalid input: weights must have the same dimensions as coefficients"
                        .to_string(),
                ));
            }
        }
        if num_rows <= degree_u {
            return Err(AlgebraError::new(format!(
                "NurbsSurface invalid input: number of rows ({}) <= degree_u ({})",
                num_rows, degree_u
            )));
        }
        if num_cols <= degree_v {
            return Err(AlgebraError::new(format!(
                "NurbsSurface invalid input: number of columns ({}) <= degree_v ({})",
                num_cols, degree_v
            )));
        }
        if knot_vector_u.len() != num_rows + 1 + degree_u {
            return Err(AlgebraError::new(format!(
                "NurbsSurface invalid input: knot_vector_u.len() ({}) != number of rows ({}) + 1 + degree_u ({})",
                knot_vector_u.len(),
                num_rows,
                degree_u
            )));
        }
        if knot_vector_v.len() != num_cols + 1 + degree_v {
            return Err(AlgebraError::new(format!(
                "NurbsSurface invalid input: knot_vector_v.len() ({}) != number of columns ({}) + 1 + degree_v ({})",
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
            weights,
            knot_vector_u,
            knot_vector_v,
            degree_u,
            degree_v,
        })
    }

    /// Constructs a NurbsSurface corresponding to a tensor–product basis function.
    ///
    /// All control points are set to zero except for the one at (index_u, index_v) which is set to
    /// `unit_vector` and its weight set to one.
    pub fn try_new_from_basis(
        index_u: usize,
        index_v: usize,
        degree_u: usize,
        degree_v: usize,
        knot_vector_u: Vec<EFloat64>,
        knot_vector_v: Vec<EFloat64>,
        unit_vector: Point,
    ) -> AlgebraResult<NurbsSurface> {
        let num_rows = knot_vector_u.len() - 1 - degree_u;
        let num_cols = knot_vector_v.len() - 1 - degree_v;
        if index_u >= num_rows {
            return Err(AlgebraError::new(format!(
                "NurbsSurface invalid input: index_u {} is out of bounds (num_rows = {})",
                index_u, num_rows
            )));
        }
        if index_v >= num_cols {
            return Err(AlgebraError::new(format!(
                "NurbsSurface invalid input: index_v {} is out of bounds (num_cols = {})",
                index_v, num_cols
            )));
        }
        let mut coefficients = vec![vec![Point::zero(); num_cols]; num_rows];
        let mut weights = vec![vec![EFloat64::zero(); num_cols]; num_rows];
        coefficients[index_u][index_v] = unit_vector;
        weights[index_u][index_v] = EFloat64::one();
        Self::try_new(
            coefficients,
            weights,
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

    /// Finds the knot span index for a given parameter `t` in the u–direction.
    pub fn find_span_u(&self, t: EFloat64) -> Option<usize> {
        NurbsSurface::find_span_generic(&self.knot_vector_u, t)
    }

    /// Finds the knot span index for a given parameter `t` in the v–direction.
    pub fn find_span_v(&self, t: EFloat64) -> Option<usize> {
        NurbsSurface::find_span_generic(&self.knot_vector_v, t)
    }

    /// Evaluates the NURBS surface at the parameter pair (u, v).
    ///
    /// This is done by first converting the control net into homogeneous coordinates (multiplying
    /// each control point by its weight and augmenting with the weight), then applying the two–step
    /// de Boor algorithm (first along the v–direction for each row, then along the u–direction),
    /// and finally converting back to Euclidean space.
    pub fn eval(&self, u: EFloat64, v: EFloat64) -> Point {
        let k_u = match self.find_span_u(u.clone()) {
            Some(span) => span,
            None => return Point::zero(),
        };
        let k_v = match self.find_span_v(v.clone()) {
            Some(span) => span,
            None => return Point::zero(),
        };
        let p_u = self.degree_u;
        let p_v = self.degree_v;

        // Build a 2D array of homogeneous control points.
        // Each homogeneous point is of the form (w * P, w).
        let mut d: Vec<Vec<NurbsHelperPoint>> = Vec::with_capacity(p_u + 1);
        for i in 0..=p_u {
            if k_u + i < p_u || k_u + i - p_u >= self.coefficients.len() {
                d.push(vec![NurbsHelperPoint::zero(); self.coefficients[0].len()]);
            } else {
                let mut row = Vec::with_capacity(p_v + 1);
                for j in 0..=p_v {
                    if k_v + j < p_v || k_v + j - p_v >= self.coefficients[0].len() {
                        row.push(NurbsHelperPoint::zero());
                    } else {
                        row.push(NurbsHelperPoint {
                            point: self.coefficients[k_u + i - p_u][k_v + j - p_v]
                                * self.weights[k_u + i - p_u][k_v + j - p_v],
                            weight: self.weights[k_u + i - p_u][k_v + j - p_v].clone(),
                        });
                    }
                }
                d.push(row);
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

        let dh = d[p_u][p_v].clone();
        (dh.point / dh.weight).unwrap_or(Point::zero())
    }
}

impl Display for NurbsSurface {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "NurbsSurface(")?;
        writeln!(f, "  Control Points:")?;
        for row in &self.coefficients {
            write!(f, "    [")?;
            for pt in row {
                write!(f, "{}, ", pt)?;
            }
            writeln!(f, "]")?;
        }
        writeln!(f, "  Weights:")?;
        for row in &self.weights {
            write!(f, "    [")?;
            for w in row {
                write!(f, "{}, ", w)?;
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

/// Helper struct for representing homogeneous control points during evaluation.
#[derive(Clone)]
struct NurbsHelperPoint {
    point: Point,
    weight: EFloat64,
}

impl NurbsHelperPoint {
    pub fn zero() -> Self {
        Self {
            point: Point::zero(),
            weight: EFloat64::zero(),
        }
    }
}

impl std::ops::Add for NurbsHelperPoint {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            point: self.point + other.point,
            weight: self.weight + other.weight,
        }
    }
}

impl std::ops::Mul<EFloat64> for NurbsHelperPoint {
    type Output = Self;
    fn mul(self, scalar: EFloat64) -> Self {
        Self {
            point: self.point * scalar,
            weight: self.weight * scalar,
        }
    }
}

impl std::ops::Mul<NurbsHelperPoint> for EFloat64 {
    type Output = NurbsHelperPoint;
    fn mul(self, rhs: NurbsHelperPoint) -> NurbsHelperPoint {
        rhs * self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::efloat::EFloat64;
    use crate::primitives::point::Point;

    fn to_efloat_vec(values: Vec<f64>) -> Vec<EFloat64> {
        values.into_iter().map(EFloat64::from).collect()
    }

    #[test]
    fn test_nurbs_surface_eval_unit_weights() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with all weights equal to one.
        // For clarity, each control point is set to a scalar multiple of unit_z.
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
        let weights = vec![
            vec![EFloat64::one(), EFloat64::one(), EFloat64::one()],
            vec![EFloat64::one(), EFloat64::one(), EFloat64::one()],
            vec![EFloat64::one(), EFloat64::one(), EFloat64::one()],
        ];

        // Use clamped knot vectors for a quadratic (degree 2) surface.
        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;

        // For a clamped surface, the evaluation at the boundaries should equal the corresponding control points.
        // Evaluate at (u,v) = (0,0) → expected to equal the control point at (0,0)
        let p00 = surface.eval(EFloat64::from(0.0), EFloat64::from(0.0));
        assert_eq!(p00, Point::unit_z() * EFloat64::from(1.0));

        // Evaluate at (u,v) = (1,0) → expected to equal the control point at (2,0) (last row in u–direction)
        let p10 = surface.eval(EFloat64::from(0.9999), EFloat64::from(0.0));
        assert_eq!(p10, Point::unit_z() * EFloat64::new(7.0, 6.999));

        // Evaluate at (u,v) = (0,1) → expected to equal the control point at (0,2) (last column in v–direction)
        let p01 = surface.eval(EFloat64::from(0.0), EFloat64::from(0.9999));
        assert_eq!(p01, Point::unit_z() * EFloat64::new(3.0, 2.999));

        // Evaluate at (u,v) = (1,1) → expected to equal the control point at (2,2)
        let p11 = surface.eval(EFloat64::from(0.9999), EFloat64::from(0.9999));
        assert_eq!(p11, Point::unit_z() * EFloat64::new(9.0, 8.999));

        // Evaluate at an interior parameter value.
        let pint = surface.eval(EFloat64::from(0.5), EFloat64::from(0.5));
        // We do not have a hand-calculated expected value here, so check it is nonzero.
        assert!(pint != Point::zero());
        Ok(())
    }

    #[test]
    fn test_nurbs_surface_eval_non_unit_weights() -> AlgebraResult<()> {
        // Create a 3x3 NURBS surface with varying weights.
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
        let weights = vec![
            vec![
                EFloat64::from(1.0),
                EFloat64::from(2.0),
                EFloat64::from(1.0),
            ],
            vec![
                EFloat64::from(1.0),
                EFloat64::from(3.0),
                EFloat64::from(1.0),
            ],
            vec![
                EFloat64::from(1.0),
                EFloat64::from(2.0),
                EFloat64::from(1.0),
            ],
        ];

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);

        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;

        // For a clamped surface, the evaluation at the boundaries should still equal the corresponding control points.
        let p00 = surface.eval(EFloat64::from(0.0), EFloat64::from(0.0));
        let p10 = surface.eval(EFloat64::from(0.99999), EFloat64::from(0.0));
        let p01 = surface.eval(EFloat64::from(0.0), EFloat64::from(0.99999));
        let p11 = surface.eval(EFloat64::from(0.99999), EFloat64::from(0.99999));

        // Despite non-uniform weights, the evaluation at the clamped ends equals the control points.
        assert_eq!(p00, Point::unit_z() * EFloat64::from(1.0));
        assert_eq!(p10, Point::unit_z() * EFloat64::new(7.0, 6.999));
        assert_eq!(p01, Point::unit_z() * EFloat64::new(3.0, 2.999));
        assert_eq!(p11, Point::unit_z() * EFloat64::new(9.0, 8.999));

        let pint = surface.eval(EFloat64::from(0.5), EFloat64::from(0.5));
        assert!(pint != Point::zero());
        Ok(())
    }

    #[test]
    fn test_nurbs_surface_from_basis() -> AlgebraResult<()> {
        // Test constructing a NurbsSurface using try_new_from_basis.
        // Here, we create a 3x3 surface which is zero everywhere except at the control point at (1,1).
        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface = NurbsSurface::try_new_from_basis(
            1,
            1,
            2,
            2,
            knot_vector_u,
            knot_vector_v,
            Point::unit_z(),
        )?;

        // Since only the control point at (1,1) is nonzero (with weight one),
        // evaluate the surface at several parameters.
        let p = surface.eval(EFloat64::from(0.5), EFloat64::from(0.5));
        // We expect a nonzero value because the basis function centered at (1,1) has support in the interior.
        assert!(p != Point::zero());
        Ok(())
    }
}
