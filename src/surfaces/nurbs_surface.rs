use std::fmt::Display;

use crate::{
    algebra_error::{AlgebraError, AlgebraResult},
    curves::nurbs_curve::NurbsCurve,
    primitives::nurb_helper_point::NurbHelperPoint,
    primitives::{convex_hull::ConvexHull, efloat::EFloat64, point::Point},
    surfaces::surface_like::SurfaceLike,
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

    /// Inserts the knot value `t` once into the NURBS surface in the u–direction using the standard
    /// knot insertion algorithm. This updates the u knot vector, control net, and weights.
    pub fn insert_knot_u(&self, t: EFloat64) -> AlgebraResult<NurbsSurface> {
        let k = match self.find_span_u(t.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Parameter t is out of the valid domain for knot insertion in u".to_string(),
                ));
            }
        };
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

        // Build new control net and weights (u–direction only changes: columns remain unchanged).
        let mut new_coeffs: Vec<Vec<Point>> = Vec::with_capacity(self.coefficients.len() + 1);
        let mut new_weights: Vec<Vec<EFloat64>> = Vec::with_capacity(self.weights.len() + 1);

        // Copy rows that are not affected.
        for i in 0..(k - p + 1) {
            new_coeffs.push(self.coefficients[i].clone());
            new_weights.push(self.weights[i].clone());
        }

        // Recompute affected rows.
        for i in (k - p + 1)..=k {
            let mut new_row_coeffs: Vec<Point> = Vec::with_capacity(num_cols);
            let mut new_row_weights: Vec<EFloat64> = Vec::with_capacity(num_cols);
            let alpha = ((t.clone() - self.knot_vector_u[i].clone())
                / (self.knot_vector_u[i + p].clone() - self.knot_vector_u[i].clone()))
            .unwrap_or(EFloat64::zero());
            for j in 0..num_cols {
                let weight = self.weights[i - 1][j].clone() * (EFloat64::one() - alpha.clone())
                    + self.weights[i][j].clone() * alpha.clone();
                let alpha = self.coefficients[i - 1][j].clone()
                    * self.weights[i - 1][j].clone()
                    * (EFloat64::one() - alpha.clone())
                    + self.coefficients[i][j].clone() * self.weights[i][j].clone() * alpha.clone();
                let pt = match alpha.clone() / weight.clone() {
                    Ok(pt) => pt,
                    Err(_) => Point::zero(),
                };
                new_row_coeffs.push(pt);
                new_row_weights.push(weight);
            }
            new_coeffs.push(new_row_coeffs);
            new_weights.push(new_row_weights);
        }

        // Copy the remaining rows.
        for i in k..=m {
            new_coeffs.push(self.coefficients[i].clone());
            new_weights.push(self.weights[i].clone());
        }

        NurbsSurface::try_new(
            new_coeffs,
            new_weights,
            new_knot_vector_u,
            self.knot_vector_v.clone(),
            p,
            self.degree_v,
        )
    }

    /// Inserts the knot value `t` once into the NURBS surface in the v–direction using the standard
    /// knot insertion algorithm. This updates the v knot vector, control net, and weights.
    pub fn insert_knot_v(&self, t: EFloat64) -> AlgebraResult<NurbsSurface> {
        let k = match self.find_span_v(t.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Parameter t is out of the valid domain for knot insertion in v".to_string(),
                ));
            }
        };
        let p = self.degree_v;
        let n = self.coefficients[0].len() - 1; // last column index
        let num_rows = self.coefficients.len();

        // Build new v knot vector.
        let mut new_knot_vector_v: Vec<EFloat64> = Vec::with_capacity(self.knot_vector_v.len() + 1);
        for i in 0..=k {
            new_knot_vector_v.push(self.knot_vector_v[i].clone());
        }
        new_knot_vector_v.push(t.clone());
        for i in (k + 1)..self.knot_vector_v.len() {
            new_knot_vector_v.push(self.knot_vector_v[i].clone());
        }

        // Build new control net and weights (v–direction only changes: rows remain unchanged).
        let mut new_coeffs: Vec<Vec<Point>> = Vec::with_capacity(num_rows);
        let mut new_weights: Vec<Vec<EFloat64>> = Vec::with_capacity(num_rows);
        for i in 0..num_rows {
            let mut new_row_coeffs: Vec<Point> = Vec::with_capacity(self.coefficients[i].len() + 1);
            let mut new_row_weights: Vec<EFloat64> = Vec::with_capacity(self.weights[i].len() + 1);
            // Copy columns that are not affected.
            for j in 0..(k - p + 1) {
                new_row_coeffs.push(self.coefficients[i][j].clone());
                new_row_weights.push(self.weights[i][j].clone());
            }
            // Recompute affected columns.
            for j in (k - p + 1)..=k {
                let alpha = ((t.clone() - self.knot_vector_v[j].clone())
                    / (self.knot_vector_v[j + p].clone() - self.knot_vector_v[j].clone()))
                .unwrap_or(EFloat64::zero());
                let weight = self.weights[i][j - 1].clone() * (EFloat64::one() - alpha.clone())
                    + self.weights[i][j].clone() * alpha.clone();
                let alpha = self.coefficients[i][j - 1].clone()
                    * self.weights[i][j - 1].clone()
                    * (EFloat64::one() - alpha.clone())
                    + self.coefficients[i][j].clone() * self.weights[i][j].clone() * alpha.clone();
                let pt = match alpha.clone() / weight.clone() {
                    Ok(pt) => pt,
                    Err(_) => Point::zero(),
                };
                new_row_coeffs.push(pt);
                new_row_weights.push(weight);
            }
            // Copy the remaining columns.
            for j in k..=n {
                new_row_coeffs.push(self.coefficients[i][j].clone());
                new_row_weights.push(self.weights[i][j].clone());
            }
            new_coeffs.push(new_row_coeffs);
            new_weights.push(new_row_weights);
        }

        NurbsSurface::try_new(
            new_coeffs,
            new_weights,
            self.knot_vector_u.clone(),
            new_knot_vector_v,
            self.degree_u,
            p,
        )
    }

    /// Subdivides the NURBS surface along the u–direction at parameter `t` into two new
    /// NurbsSurface segments. The method inserts `t` repeatedly until its multiplicity equals
    /// degree_u+1 (i.e. a break point) then splits the control net, weights, and u knot vector.
    ///
    /// The v–direction remains unchanged.
    fn subdivide_u_impl(&self, t: EFloat64) -> AlgebraResult<(NurbsSurface, NurbsSurface)> {
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
        let left_weights = surface.weights[..=(t_index - p)].to_vec();
        let left_knots_u = surface.knot_vector_u[..=t_index + 1].to_vec();
        let left_surface = NurbsSurface::try_new(
            left_coeffs,
            left_weights,
            left_knots_u,
            surface.knot_vector_v.clone(),
            p,
            surface.degree_v,
        )?;

        // Build right surface: rows from (t_index - p) to end.
        let right_coeffs = surface.coefficients[(t_index - p)..].to_vec();
        let right_weights = surface.weights[(t_index - p)..].to_vec();
        let right_knots_u = surface.knot_vector_u[t_index - p..].to_vec();
        let right_surface = NurbsSurface::try_new(
            right_coeffs,
            right_weights,
            right_knots_u,
            surface.knot_vector_v.clone(),
            p,
            surface.degree_v,
        )?;

        Ok((left_surface, right_surface))
    }

    /// Subdivides the NURBS surface along the v–direction at parameter `t` into two new
    /// NurbsSurface segments. The method inserts `t` repeatedly until its multiplicity equals
    /// degree_v+1 (i.e. a break point) then splits the control net, weights, and v knot vector.
    ///
    /// The u–direction remains unchanged.
    fn subdivide_v_impl(&self, t: EFloat64) -> AlgebraResult<(NurbsSurface, NurbsSurface)> {
        let p = self.degree_v;

        if t < self.knot_vector_v[p] || t > self.knot_vector_v[self.knot_vector_v.len() - p - 1] {
            return Err(AlgebraError::new(
                "Parameter t is out of the valid domain for v subdivision".to_string(),
            ));
        }

        // Determine current multiplicity of t in v–direction.
        let current_multiplicity = self.knot_vector_v.iter().filter(|&knot| *knot == t).count();
        // t must appear with multiplicity p+1.
        let r = p - current_multiplicity;
        let mut surface = self.clone();
        for _ in 0..r {
            surface = surface.insert_knot_v(t.clone())?;
        }

        let t_index = match surface.find_span_v(t.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Failed to locate knot t with full multiplicity after insertion in v"
                        .to_string(),
                ));
            }
        };

        // Build left surface: columns 0 to (t_index - p) and the corresponding v–knot vector.
        let mut left_coeffs: Vec<Vec<Point>> = Vec::with_capacity(surface.coefficients.len());
        let mut left_weights: Vec<Vec<EFloat64>> = Vec::with_capacity(surface.weights.len());
        for (i, row) in surface.coefficients.iter().enumerate() {
            left_coeffs.push(row[..=(t_index - p)].to_vec());
            left_weights.push(surface.weights[i][..=(t_index - p)].to_vec());
        }
        let left_knots_v = surface.knot_vector_v[..=t_index + 1].to_vec();
        let left_surface = NurbsSurface::try_new(
            left_coeffs,
            left_weights,
            surface.knot_vector_u.clone(),
            left_knots_v,
            surface.degree_u,
            p,
        )?;

        // Build right surface: columns from (t_index - p) to end.
        let mut right_coeffs: Vec<Vec<Point>> = Vec::with_capacity(surface.coefficients.len());
        let mut right_weights: Vec<Vec<EFloat64>> = Vec::with_capacity(surface.weights.len());
        for (i, row) in surface.coefficients.iter().enumerate() {
            right_coeffs.push(row[(t_index - p)..].to_vec());
            right_weights.push(surface.weights[i][(t_index - p)..].to_vec());
        }
        let right_knots_v = surface.knot_vector_v[t_index - p..].to_vec();
        let right_surface = NurbsSurface::try_new(
            right_coeffs,
            right_weights,
            surface.knot_vector_u.clone(),
            right_knots_v,
            surface.degree_u,
            p,
        )?;

        Ok((left_surface, right_surface))
    }

    /// Evaluates the NURBS surface at the parameter pair (u, v).
    ///
    /// This is done by first converting the control net into homogeneous coordinates (multiplying
    /// each control point by its weight and augmenting with the weight), then applying the two–step
    /// de Boor algorithm (first along the v–direction for each row, then along the u–direction),
    /// and finally converting back to Euclidean space.
    pub fn eval_impl(&self, u: EFloat64, v: EFloat64) -> Point {
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
        let mut d: Vec<Vec<NurbHelperPoint>> = Vec::with_capacity(p_u + 1);
        for i in 0..=p_u {
            if k_u + i < p_u || k_u + i - p_u >= self.coefficients.len() {
                d.push(vec![NurbHelperPoint::zero(); self.coefficients[0].len()]);
            } else {
                let mut row = Vec::with_capacity(p_v + 1);
                for j in 0..=p_v {
                    if k_v + j < p_v || k_v + j - p_v >= self.coefficients[0].len() {
                        row.push(NurbHelperPoint::zero());
                    } else {
                        row.push(NurbHelperPoint::new(
                            self.coefficients[k_u + i - p_u][k_v + j - p_v]
                                * self.weights[k_u + i - p_u][k_v + j - p_v],
                            self.weights[k_u + i - p_u][k_v + j - p_v].clone(),
                        ));
                    }
                }
                d.push(row);
            }
        }

        // Apply de Boor's algorithm.
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
        dh.to_point().unwrap_or(Point::zero())
    }

    /// Returns a NURBS curve that represents the iso-curve of the surface at a fixed u parameter.
    /// The resulting curve is defined over the v parameter domain.
    pub fn iso_curve_u(&self, u: EFloat64) -> AlgebraResult<NurbsCurve> {
        let k_u = match self.find_span_u(u.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Parameter u is out of the valid domain for iso-curve extraction".to_string(),
                ));
            }
        };
        let p_u = self.degree_u;

        // Initialize arrays for the de Boor algorithm
        let mut d: Vec<Vec<NurbHelperPoint>> = Vec::with_capacity(self.coefficients[0].len());
        for j in 0..self.coefficients[0].len() {
            let mut row = Vec::with_capacity(p_u + 1);
            for i in 0..=p_u {
                if k_u + i < p_u || k_u + i - p_u >= self.coefficients.len() {
                    row.push(NurbHelperPoint::zero());
                } else {
                    row.push(NurbHelperPoint::new(
                        self.coefficients[k_u + i - p_u][j].clone()
                            * self.weights[k_u + i - p_u][j].clone(),
                        self.weights[k_u + i - p_u][j].clone(),
                    ));
                }
            }
            d.push(row);
        }

        // Apply de Boor's algorithm in the u-direction
        for r_u in 1..=p_u {
            for j_u in (r_u..=p_u).rev() {
                let alpha_u =
                    match k_u + j_u < p_u || j_u + 1 + k_u - r_u >= self.knot_vector_u.len() {
                        true => EFloat64::zero(),
                        false => {
                            let left_knot = self.knot_vector_u[j_u + k_u - p_u].clone();
                            let right_knot = self.knot_vector_u[j_u + 1 + k_u - r_u].clone();
                            ((u.clone() - left_knot) / (right_knot - left_knot))
                                .unwrap_or(EFloat64::zero())
                        }
                    };
                for j in 0..self.coefficients[0].len() {
                    d[j][j_u] = d[j][j_u - 1].clone() * (EFloat64::one() - alpha_u.clone())
                        + d[j][j_u].clone() * alpha_u.clone();
                }
            }
        }

        // Extract control points and weights for the resulting NURBS curve
        let mut curve_coeffs = Vec::with_capacity(self.coefficients[0].len());
        let mut curve_weights = Vec::with_capacity(self.coefficients[0].len());
        for j in 0..self.coefficients[0].len() {
            let dh = d[j][p_u].clone();
            match dh.to_point() {
                Ok(pt) => curve_coeffs.push(pt),
                Err(_) => curve_coeffs.push(Point::zero()),
            }
            curve_weights.push(d[j][p_u].weight());
        }

        NurbsCurve::try_new(
            curve_coeffs,
            curve_weights,
            self.knot_vector_v.clone(),
            self.degree_v,
        )
    }

    /// Returns a NURBS curve that represents the iso-curve of the surface at a fixed v parameter.
    /// The resulting curve is defined over the u parameter domain.
    pub fn iso_curve_v(&self, v: EFloat64) -> AlgebraResult<NurbsCurve> {
        let k_v = match self.find_span_v(v.clone()) {
            Some(idx) => idx,
            None => {
                return Err(AlgebraError::new(
                    "Parameter v is out of the valid domain for iso-curve extraction".to_string(),
                ));
            }
        };
        let p_v = self.degree_v;

        // Initialize arrays for the de Boor algorithm
        let mut d: Vec<Vec<NurbHelperPoint>> = Vec::with_capacity(self.coefficients.len());
        for i in 0..self.coefficients.len() {
            let mut row = Vec::with_capacity(p_v + 1);
            for j in 0..=p_v {
                if k_v + j < p_v || k_v + j - p_v >= self.coefficients[0].len() {
                    row.push(NurbHelperPoint::zero());
                } else {
                    row.push(NurbHelperPoint::new(
                        self.coefficients[i][k_v + j - p_v].clone()
                            * self.weights[i][k_v + j - p_v].clone(),
                        self.weights[i][k_v + j - p_v].clone(),
                    ));
                }
            }
            d.push(row);
        }

        // Apply de Boor's algorithm in the v-direction
        for r_v in 1..=p_v {
            for j_v in (r_v..=p_v).rev() {
                let alpha_v =
                    match k_v + j_v < p_v || j_v + 1 + k_v - r_v >= self.knot_vector_v.len() {
                        true => EFloat64::zero(),
                        false => {
                            let left_knot = self.knot_vector_v[j_v + k_v - p_v].clone();
                            let right_knot = self.knot_vector_v[j_v + 1 + k_v - r_v].clone();
                            ((v.clone() - left_knot) / (right_knot - left_knot))
                                .unwrap_or(EFloat64::zero())
                        }
                    };
                for i in 0..self.coefficients.len() {
                    d[i][j_v] = d[i][j_v - 1].clone() * (EFloat64::one() - alpha_v.clone())
                        + d[i][j_v].clone() * alpha_v.clone();
                }
            }
        }

        // Extract control points and weights for the resulting NURBS curve
        let mut curve_coeffs = Vec::with_capacity(self.coefficients.len());
        let mut curve_weights = Vec::with_capacity(self.coefficients.len());
        for i in 0..self.coefficients.len() {
            let dh = d[i][p_v].clone();
            match dh.to_point() {
                Ok(pt) => curve_coeffs.push(pt),
                Err(_) => curve_coeffs.push(Point::zero()),
            }
            curve_weights.push(d[i][p_v].weight());
        }

        NurbsCurve::try_new(
            curve_coeffs,
            curve_weights,
            self.knot_vector_u.clone(),
            self.degree_u,
        )
    }

    /// Returns the convex hull of the control points of the NURBS surface.
    /// This is the smallest convex set containing all control points.
    pub fn get_convex_hull(&self) -> AlgebraResult<ConvexHull> {
        // Flatten the 2D control points array into a 1D vector
        let points: Vec<Point> = self
            .coefficients
            .iter()
            .flat_map(|row| row.iter().cloned())
            .collect();

        // Use the existing ConvexHull implementation
        ConvexHull::try_new(points)
    }
}

impl SurfaceLike for NurbsSurface {
    fn u_span(&self) -> (EFloat64, EFloat64) {
        (
            self.knot_vector_u[self.degree_u].clone(),
            self.knot_vector_u[self.knot_vector_u.len() - self.degree_u - 1].clone(),
        )
    }

    fn v_span(&self) -> (EFloat64, EFloat64) {
        (
            self.knot_vector_v[self.degree_v].clone(),
            self.knot_vector_v[self.knot_vector_v.len() - self.degree_v - 1].clone(),
        )
    }

    fn eval(&self, u: EFloat64, v: EFloat64) -> Point {
        self.eval_impl(u, v)
    }

    fn subdivide_u(
        &self,
        t: EFloat64,
    ) -> AlgebraResult<(Box<dyn SurfaceLike>, Box<dyn SurfaceLike>)> {
        let (left, right) = self.subdivide_u_impl(t)?;
        Ok((Box::new(left), Box::new(right)))
    }

    fn subdivide_v(
        &self,
        t: EFloat64,
    ) -> AlgebraResult<(Box<dyn SurfaceLike>, Box<dyn SurfaceLike>)> {
        let (left, right) = self.subdivide_v_impl(t)?;
        Ok((Box::new(left), Box::new(right)))
    }

    fn get_convex_hull(&self) -> AlgebraResult<ConvexHull> {
        // Flatten the 2D control points array into a 1D vector
        let points: Vec<Point> = self
            .coefficients
            .iter()
            .flat_map(|row| row.iter().cloned())
            .collect();

        // Use the existing ConvexHull implementation
        ConvexHull::try_new(points)
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

    #[test]
    fn test_insert_knot_u() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with varying weights.
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

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;
        println!("Original surface: {}", surface);

        // Insert a knot at u = 1.0
        let t = EFloat64::from(1.0);
        let surface2 = surface.insert_knot_u(t.clone())?;
        println!("Surface after inserting knot u at {}: {}", t, surface2);

        // Check that evaluation remains the same after knot insertion.
        // We test at various points in the parameter domain.
        for i in 0..=100 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            for j in 0..=100 {
                let v_val = EFloat64::from(j as f64 / 100.0);
                let orig = surface.eval(u_val.clone(), v_val.clone());
                let new = surface2.eval(u_val.clone(), v_val.clone());
                assert_eq!(orig, new, "Mismatch at (u,v)=({}, {})", u_val, v_val);
            }
        }

        // Verify that the knot vector was updated correctly
        assert_eq!(
            surface2.knot_vector_u.len(),
            surface.knot_vector_u.len() + 1
        );
        assert!(surface2.knot_vector_u.contains(&t));

        // Verify that the control net and weights were updated correctly
        assert_eq!(surface2.coefficients.len(), surface.coefficients.len() + 1);
        assert_eq!(surface2.weights.len(), surface.weights.len() + 1);

        Ok(())
    }

    #[test]
    fn test_subdivide_u() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with varying weights.
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

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;

        // Subdivide at u = 1.0
        let t = EFloat64::from(1.0);
        let (left, right) = surface.subdivide_u(t.clone())?;

        println!("Left: {}", left);
        println!("Right: {}", right);

        // Verify that evaluation remains the same after subdivision
        // Test points in the left segment (u < 1.0)
        for i in 0..=50 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            if u_val < t {
                for j in 0..=100 {
                    let v_val = EFloat64::from(j as f64 / 100.0);
                    let orig = surface.eval(u_val.clone(), v_val.clone());
                    let new = left.eval(u_val.clone(), v_val.clone());
                    assert_eq!(
                        orig, new,
                        "Mismatch in left segment at (u,v)=({}, {})",
                        u_val, v_val
                    );
                }
            }
        }

        // Test points in the right segment (u > 1.0)
        for i in 51..=100 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            if u_val > t {
                for j in 0..=100 {
                    let v_val = EFloat64::from(j as f64 / 100.0);
                    let orig = surface.eval(u_val.clone(), v_val.clone());
                    let new = right.eval(u_val.clone(), v_val.clone());
                    assert_eq!(
                        orig, new,
                        "Mismatch in right segment at (u,v)=({}, {})",
                        u_val, v_val
                    );
                }
            }
        }

        // Verify that the knot vectors were updated correctly
        let (left_min_u, left_max_u) = left.u_span();
        let (right_min_u, right_max_u) = right.u_span();
        assert!(left_min_u <= t && t <= left_max_u);
        assert!(right_min_u <= t && t <= right_max_u);

        Ok(())
    }

    #[test]
    fn test_insert_knot_v() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with varying weights.
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

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]); // Changed to [0,2] domain
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;
        println!("Original surface: {}", surface);

        // Insert a knot at v = 1.0 (now in the middle of the domain)
        let t = EFloat64::from(1.0);
        let surface2 = surface.insert_knot_v(t.clone())?;
        println!("Surface after inserting knot v at {}: {}", t, surface2);

        // Check that evaluation remains the same after knot insertion.
        // We test at various points in the parameter domain.
        for i in 0..=100 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            for j in 0..=100 {
                let v_val = EFloat64::from(j as f64 / 100.0 * 2.0); // Scale to [0,2] domain
                let orig = surface.eval(u_val.clone(), v_val.clone());
                let new = surface2.eval(u_val.clone(), v_val.clone());
                assert_eq!(orig, new, "Mismatch at (u,v)=({}, {})", u_val, v_val);
            }
        }

        // Verify that the knot vector was updated correctly
        assert_eq!(
            surface2.knot_vector_v.len(),
            surface.knot_vector_v.len() + 1
        );
        assert!(surface2.knot_vector_v.contains(&t));

        // Verify that the control net and weights were updated correctly
        assert_eq!(surface2.coefficients.len(), surface.coefficients.len());
        assert_eq!(surface2.weights.len(), surface.weights.len());

        Ok(())
    }

    #[test]
    fn test_subdivide_v() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with varying weights
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

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]); // Changed to [0,2] domain
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;

        // Subdivide at v = 1.0 (now in the middle of the domain)
        let t = EFloat64::from(1.0);
        let (left, right) = surface.subdivide_v(t.clone())?;

        println!("Left: {}", left);
        println!("Right: {}", right);

        // Verify that evaluation remains the same after subdivision
        // Test points in the left segment (v < 1.0)
        for i in 0..=100 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            for j in 0..=50 {
                let v_val = EFloat64::from(j as f64 / 100.0 * 2.0);
                if v_val < t {
                    let orig = surface.eval(u_val.clone(), v_val.clone());
                    let new = left.eval(u_val.clone(), v_val.clone());
                    assert_eq!(
                        orig, new,
                        "Mismatch in left segment at (u,v)=({}, {})",
                        u_val, v_val
                    );
                }
            }
        }

        // Test points in the right segment (v > 1.0)
        for i in 0..=100 {
            let u_val = EFloat64::from(i as f64 / 100.0 * 2.0);
            for j in 51..=100 {
                let v_val = EFloat64::from(j as f64 / 100.0 * 2.0);
                if v_val > t {
                    let orig = surface.eval(u_val.clone(), v_val.clone());
                    let new = right.eval(u_val.clone(), v_val.clone());
                    assert_eq!(
                        orig, new,
                        "Mismatch in right segment at (u,v)=({}, {})",
                        u_val, v_val
                    );
                }
            }
        }

        // Verify that the knot vectors were updated correctly
        let (left_min_v, left_max_v) = left.v_span();
        let (right_min_v, right_max_v) = right.v_span();
        assert!(left_min_v <= t && t <= left_max_v);
        assert!(right_min_v <= t && t <= right_max_v);

        Ok(())
    }

    #[test]
    fn test_iso_curves() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with varying weights
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

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0]);
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;

        // Test u-iso-curve at u = 1.0
        let u_fixed = EFloat64::from(0.4523453);
        let u_iso = surface.iso_curve_u(u_fixed.clone())?;
        // Verify that evaluating the iso-curve at any v gives the same result as evaluating the surface
        for i in 0..=100 {
            let v = EFloat64::from(i as f64 / 100.0 * 2.0);
            let surface_point = surface.eval(u_fixed.clone(), v.clone());
            let curve_point = u_iso.eval(v);
            assert_eq!(
                surface_point, curve_point,
                "Mismatch in u-iso-curve at v={}",
                v
            );
        }

        // Test v-iso-curve at v = 1.0
        let v_fixed = EFloat64::from(0.4523453);
        let v_iso = surface.iso_curve_v(v_fixed.clone())?;
        // Verify that evaluating the iso-curve at any u gives the same result as evaluating the surface
        for i in 0..=100 {
            let u = EFloat64::from(i as f64 / 100.0 * 2.0);
            let surface_point = surface.eval(u.clone(), v_fixed.clone());
            let curve_point = v_iso.eval(u);
            assert_eq!(
                surface_point, curve_point,
                "Mismatch in v-iso-curve at u={}",
                u
            );
        }

        // Test error cases
        // Try to get iso-curve outside the valid domain
        assert!(surface.iso_curve_u(EFloat64::from(-1.0)).is_err());
        assert!(surface.iso_curve_u(EFloat64::from(3.0)).is_err());
        assert!(surface.iso_curve_v(EFloat64::from(-1.0)).is_err());
        assert!(surface.iso_curve_v(EFloat64::from(3.0)).is_err());

        Ok(())
    }

    #[test]
    fn test_convex_hull_contains_all_points() -> AlgebraResult<()> {
        // Create a simple 3x3 NURBS surface with control points forming a pyramid-like shape
        let coefficients = vec![
            vec![
                Point::new(
                    EFloat64::from(0.0),
                    EFloat64::from(0.0),
                    EFloat64::from(0.0),
                ),
                Point::new(
                    EFloat64::from(1.0),
                    EFloat64::from(0.0),
                    EFloat64::from(0.0),
                ),
                Point::new(
                    EFloat64::from(2.0),
                    EFloat64::from(0.0),
                    EFloat64::from(0.0),
                ),
            ],
            vec![
                Point::new(
                    EFloat64::from(0.0),
                    EFloat64::from(1.0),
                    EFloat64::from(0.0),
                ),
                Point::new(
                    EFloat64::from(1.0),
                    EFloat64::from(1.0),
                    EFloat64::from(2.0),
                ), // Peak
                Point::new(
                    EFloat64::from(2.0),
                    EFloat64::from(1.0),
                    EFloat64::from(0.0),
                ),
            ],
            vec![
                Point::new(
                    EFloat64::from(0.0),
                    EFloat64::from(2.0),
                    EFloat64::from(0.0),
                ),
                Point::new(
                    EFloat64::from(1.0),
                    EFloat64::from(2.0),
                    EFloat64::from(0.0),
                ),
                Point::new(
                    EFloat64::from(2.0),
                    EFloat64::from(2.0),
                    EFloat64::from(0.0),
                ),
            ],
        ];
        let weights = vec![
            vec![EFloat64::one(), EFloat64::one(), EFloat64::one()],
            vec![EFloat64::one(), EFloat64::one(), EFloat64::one()],
            vec![EFloat64::one(), EFloat64::one(), EFloat64::one()],
        ];

        let knot_vector_u = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let knot_vector_v = to_efloat_vec(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let surface =
            NurbsSurface::try_new(coefficients, weights, knot_vector_u, knot_vector_v, 2, 2)?;

        // Get the convex hull
        let hull = surface.get_convex_hull()?;

        // Verify that all control points are inside the convex hull
        for row in &surface.coefficients {
            for point in row {
                assert!(
                    hull.contains_point(point),
                    "Control point {} is not inside the convex hull",
                    point
                );
            }
        }

        Ok(())
    }
}
