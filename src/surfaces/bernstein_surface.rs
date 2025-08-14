use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

use crate::algebra_error::AlgebraResult;
use crate::bernstein::bernstein_polynomial::BernsteinPolynomial;
use crate::binomial_coefficient;
use crate::primitives::{convex_hull::ConvexHull, efloat::EFloat64, point::Point};
use crate::surfaces::surface_like::SurfaceLike;
use crate::zero::Zero;

/// Tensor-product Bernstein surface with coefficients laid out as a 2D grid.
/// Row count is degree_u+1, column count is degree_v+1.
#[derive(Debug, Clone, PartialEq)]
pub struct BernsteinSurface<T> {
    pub coefficients: Vec<Vec<T>>, // [i_u][i_v]
}

impl<T> BernsteinSurface<T> {
    pub fn new(coefficients: Vec<Vec<T>>) -> Self {
        Self { coefficients }
    }

    pub fn degree_u(&self) -> usize {
        self.coefficients.len().saturating_sub(1)
    }

    pub fn degree_v(&self) -> usize {
        if self.coefficients.is_empty() {
            0
        } else {
            self.coefficients[0].len().saturating_sub(1)
        }
    }
}

impl<T> Display for BernsteinSurface<T>
where
    T: Display + Zero + PartialEq,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let du = self.degree_u();
        let dv = self.degree_v();
        write!(f, "BernsteinSurface(\n")?;
        for i in 0..=du {
            write!(f, "[")?;
            for j in 0..=dv {
                let c = &self.coefficients[i][j];
                // if *c != T::zero() {
                // write!(f, "{} B_{{{},{}}}(u) B_{{{},{}}}(v)", c, i, du, j, dv)?;
                write!(f, "{}", c)?;
                // }
                if j < dv {
                    write!(f, ", ")?;
                }
            }
            write!(f, "]")?;
            if i < du {
                write!(f, "\n")?;
            }
        }
        write!(f, ")")
    }
}

impl<T> BernsteinSurface<T>
where
    T: Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    /// Evaluate via tensor de Casteljau: first in u for each column, then in v; or vice versa.
    pub fn eval(&self, u: EFloat64, v: EFloat64) -> T {
        if self.coefficients.is_empty() || self.coefficients[0].is_empty() {
            panic!("Empty BernsteinSurface");
        }
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();

        // de Casteljau in u for each fixed v-index -> get a column vector along u reduced to 1 value per v
        let mut column_vals: Vec<T> = Vec::with_capacity(nv);
        for j in 0..nv {
            let mut beta: Vec<T> = (0..nu).map(|i| self.coefficients[i][j].clone()).collect();
            for r in 1..nu {
                for i in 0..(nu - r) {
                    beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                        + beta[i + 1].clone() * u.clone();
                }
            }
            column_vals.push(beta[0].clone());
        }

        // de Casteljau in v on the reduced column
        let mut gamma = column_vals;
        for r in 1..nv {
            for j in 0..(nv - r) {
                gamma[j] = gamma[j].clone() * (EFloat64::one() - v.clone())
                    + gamma[j + 1].clone() * v.clone();
            }
        }
        gamma[0].clone()
    }

    /// Extract the isoparametric curve in u for a fixed v parameter.
    /// Returns a 1D Bernstein polynomial whose control points are obtained by
    /// de Casteljau reduction along the v-direction for each u-index.
    pub fn iso_u_at(&self, v: EFloat64) -> BernsteinPolynomial<T> {
        if self.coefficients.is_empty() || self.coefficients[0].is_empty() {
            panic!("Empty BernsteinSurface");
        }
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();

        let mut control_points: Vec<T> = Vec::with_capacity(nu);
        for i in 0..nu {
            // Reduce row i along v
            let mut beta: Vec<T> = (0..nv).map(|j| self.coefficients[i][j].clone()).collect();
            for r in 1..nv {
                for j in 0..(nv - r) {
                    beta[j] = beta[j].clone() * (EFloat64::one() - v.clone())
                        + beta[j + 1].clone() * v.clone();
                }
            }
            control_points.push(beta[0].clone());
        }
        BernsteinPolynomial::new(control_points)
    }

    /// Extract the isoparametric curve in v for a fixed u parameter.
    /// Returns a 1D Bernstein polynomial whose control points are obtained by
    /// de Casteljau reduction along the u-direction for each v-index.
    pub fn iso_v_at(&self, u: EFloat64) -> BernsteinPolynomial<T> {
        if self.coefficients.is_empty() || self.coefficients[0].is_empty() {
            panic!("Empty BernsteinSurface");
        }
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();

        let mut control_points: Vec<T> = Vec::with_capacity(nv);
        for j in 0..nv {
            // Reduce column j along u
            let mut beta: Vec<T> = (0..nu).map(|i| self.coefficients[i][j].clone()).collect();
            for r in 1..nu {
                for i in 0..(nu - r) {
                    beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                        + beta[i + 1].clone() * u.clone();
                }
            }
            control_points.push(beta[0].clone());
        }
        BernsteinPolynomial::new(control_points)
    }
}

impl<T> BernsteinSurface<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    /// Degree elevation by r_u and r_v in u and v directions respectively
    pub fn elevate_degree(&self, r_u: usize, r_v: usize) -> Self {
        let du = self.degree_u();
        let dv = self.degree_v();

        // Elevate in u
        let mut elev_u: Vec<Vec<T>> = vec![vec![T::zero(); dv + 1]; du + r_u + 1];
        for i_new in 0..=du + r_u {
            for i_old in i_new.saturating_sub(r_u)..=du.min(i_new) {
                let alpha = (EFloat64::from(
                    (binomial_coefficient(r_u, i_new - i_old) * binomial_coefficient(du, i_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(du + r_u, i_new) as f64))
                .unwrap();
                for j in 0..=dv {
                    elev_u[i_new][j] =
                        elev_u[i_new][j].clone() + self.coefficients[i_old][j].clone() * alpha;
                }
            }
        }

        // Elevate in v
        let mut elev_v: Vec<Vec<T>> = vec![vec![T::zero(); dv + r_v + 1]; du + r_u + 1];
        for j_new in 0..=dv + r_v {
            for j_old in j_new.saturating_sub(r_v)..=dv.min(j_new) {
                let beta = (EFloat64::from(
                    (binomial_coefficient(r_v, j_new - j_old) * binomial_coefficient(dv, j_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(dv + r_v, j_new) as f64))
                .unwrap();
                for i in 0..=du + r_u {
                    elev_v[i][j_new] = elev_v[i][j_new].clone() + elev_u[i][j_old].clone() * beta;
                }
            }
        }

        Self::new(elev_v)
    }

    /// Subdivide in u using tensor de Casteljau at parameter t.
    pub fn subdivide_u(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();

        // Perform de Casteljau along u for each v column and record left/right polygons
        let mut left: Vec<Vec<T>> = vec![vec![T::zero(); dv]; du];
        let mut right: Vec<Vec<T>> = vec![vec![T::zero(); dv]; du];

        for j in 0..dv {
            // copy column j
            let mut beta: Vec<T> = (0..du).map(|i| self.coefficients[i][j].clone()).collect();
            left[0][j] = beta[0].clone();
            right[du - 1][j] = beta[du - 1].clone();
            for r in 1..du {
                for i in 0..(du - r) {
                    beta[i] = beta[i].clone() * (EFloat64::one() - t.clone())
                        + beta[i + 1].clone() * t.clone();
                }
                left[r][j] = beta[0].clone();
                right[du - r - 1][j] = beta[du - r - 1].clone();
            }
        }

        (Self::new(left), Self::new(right))
    }

    /// Subdivide in v using tensor de Casteljau at parameter t.
    pub fn subdivide_v(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();

        let mut left: Vec<Vec<T>> = vec![vec![T::zero(); dv]; du];
        let mut right: Vec<Vec<T>> = vec![vec![T::zero(); dv]; du];

        for i in 0..du {
            // copy row i
            let mut beta: Vec<T> = (0..dv).map(|j| self.coefficients[i][j].clone()).collect();
            left[i][0] = beta[0].clone();
            right[i][dv - 1] = beta[dv - 1].clone();
            for r in 1..dv {
                for j in 0..(dv - r) {
                    beta[j] = beta[j].clone() * (EFloat64::one() - t.clone())
                        + beta[j + 1].clone() * t.clone();
                }
                left[i][r] = beta[0].clone();
                right[i][dv - r - 1] = beta[dv - r - 1].clone();
            }
        }

        (Self::new(left), Self::new(right))
    }
}

impl<T> BernsteinSurface<T>
where
    T: Zero
        + Clone
        + Sub<Output = T>
        + Mul<EFloat64, Output = T>
        + Div<EFloat64, Output = AlgebraResult<T>>
        + PartialEq,
{
    /// Reduce degree by one in u, returning error if not exactly reducible.
    fn reduce_degree_u_once(&self) -> AlgebraResult<Self> {
        let du = self.degree_u();
        if du == 0 {
            return Err("Cannot reduce degree_u of a constant surface".into());
        }
        let dv = self.degree_v();

        let mut q: Vec<Vec<T>> = vec![vec![T::zero(); dv + 1]; du];
        // boundary: q_0,* = a_0,*
        for j in 0..=dv {
            q[0][j] = self.coefficients[0][j].clone();
        }

        let n_ef = EFloat64::from(du as f64);
        for i in 1..du {
            let denom = EFloat64::from((du - i) as f64);
            for j in 0..=dv {
                // q_i,j = (du * a_i,j - i * q_{i-1,j})/(du-i)
                let numer = self.coefficients[i][j].clone() * n_ef.clone()
                    - q[i - 1][j].clone() * EFloat64::from(i as f64);
                q[i][j] = (numer / denom.clone())?;
            }
        }

        // boundary consistency: q_{du-1,*} == a_{du,*}
        for j in 0..=dv {
            if q[du - 1][j] != self.coefficients[du][j] {
                return Err("Degree reduction in u not possible exactly".into());
            }
        }
        Ok(Self::new(q))
    }

    /// Reduce degree by one in v, returning error if not exactly reducible.
    fn reduce_degree_v_once(&self) -> AlgebraResult<Self> {
        let dv = self.degree_v();
        if dv == 0 {
            return Err("Cannot reduce degree_v of a constant surface".into());
        }
        let du = self.degree_u();

        let mut q: Vec<Vec<T>> = vec![vec![T::zero(); dv]; du + 1];
        // boundary: q_*,0 = a_*,0
        for i in 0..=du {
            q[i][0] = self.coefficients[i][0].clone();
        }

        let n_ef = EFloat64::from(dv as f64);
        for j in 1..dv {
            let denom = EFloat64::from((dv - j) as f64);
            for i in 0..=du {
                let numer = self.coefficients[i][j].clone() * n_ef.clone()
                    - q[i][j - 1].clone() * EFloat64::from(j as f64);
                q[i][j] = (numer / denom.clone())?;
            }
        }

        for i in 0..=du {
            if q[i][dv - 1] != self.coefficients[i][dv] {
                return Err("Degree reduction in v not possible exactly".into());
            }
        }
        Ok(Self::new(q))
    }

    pub fn reduce_degree(&self, r_u: usize, r_v: usize) -> AlgebraResult<Self> {
        if r_u > self.degree_u() || r_v > self.degree_v() {
            return Err("Reduction degree exceeds surface degree".into());
        }
        let mut s = self.clone();
        for _ in 0..r_u {
            s = s.reduce_degree_u_once()?;
        }
        for _ in 0..r_v {
            s = s.reduce_degree_v_once()?;
        }
        Ok(s)
    }
}

impl<T> BernsteinSurface<T>
where
    T: Zero + Clone + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    /// First-order partial derivatives as Bernstein surfaces of degrees (du-1, dv) and (du, dv-1).
    pub fn derivative_u(&self) -> BernsteinSurface<T> {
        let du = self.degree_u();
        if du == 0 {
            return BernsteinSurface::new(vec![self.coefficients[0].clone()]);
        }
        let dv = self.degree_v();
        let scale = EFloat64::from(du as f64);
        let mut coeffs = vec![vec![T::zero(); dv + 1]; du];
        for i in 0..du {
            for j in 0..=dv {
                let diff = self.coefficients[i + 1][j].clone() - self.coefficients[i][j].clone();
                coeffs[i][j] = diff * scale.clone();
            }
        }
        BernsteinSurface::new(coeffs)
    }

    pub fn derivative_v(&self) -> BernsteinSurface<T> {
        let dv = self.degree_v();
        if dv == 0 {
            let mut base = Vec::with_capacity(self.coefficients.len());
            for row in &self.coefficients {
                base.push(vec![row[0].clone()]);
            }
            return BernsteinSurface::new(base);
        }
        let du = self.degree_u();
        let scale = EFloat64::from(dv as f64);
        let mut coeffs = vec![vec![T::zero(); dv]; du + 1];
        for i in 0..=du {
            for j in 0..dv {
                let diff = self.coefficients[i][j + 1].clone() - self.coefficients[i][j].clone();
                coeffs[i][j] = diff * scale.clone();
            }
        }
        BernsteinSurface::new(coeffs)
    }
}

// Add, Sub for two Bernstein surfaces (equalize degrees by elevation)
impl<T> BernsteinSurface<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    fn equalize_degree(
        mut self,
        mut rhs: BernsteinSurface<T>,
    ) -> (BernsteinSurface<T>, BernsteinSurface<T>) {
        let (du1, dv1) = (self.degree_u(), self.degree_v());
        let (du2, dv2) = (rhs.degree_u(), rhs.degree_v());
        if du1 < du2 {
            self = self.elevate_degree(du2 - du1, 0);
        }
        if dv1 < dv2 {
            self = self.elevate_degree(0, dv2 - dv1);
        }
        let (du1, dv1) = (self.degree_u(), self.degree_v());
        if du2 < du1 {
            rhs = rhs.elevate_degree(du1 - du2, 0);
        }
        if dv2 < dv1 {
            rhs = rhs.elevate_degree(0, dv1 - dv2);
        }
        (self, rhs)
    }
}

impl<T> Add for BernsteinSurface<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let du = lhs.degree_u();
        let dv = lhs.degree_v();
        let mut coeffs = vec![vec![T::zero(); dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                coeffs[i][j] = lhs.coefficients[i][j].clone() + rhs.coefficients[i][j].clone();
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Sub for BernsteinSurface<T>
where
    T: Zero + Clone + Sub<Output = T> + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let du = lhs.degree_u();
        let dv = lhs.degree_v();
        let mut coeffs = vec![vec![T::zero(); dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                coeffs[i][j] = lhs.coefficients[i][j].clone() - rhs.coefficients[i][j].clone();
            }
        }
        Self::new(coeffs)
    }
}

// Scalar multiplication by EFloat64
impl<T> Mul<EFloat64> for BernsteinSurface<T>
where
    T: Clone + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: EFloat64) -> Self::Output {
        let mut coeffs = self.coefficients.clone();
        for row in &mut coeffs {
            for c in row {
                *c = c.clone() * rhs;
            }
        }
        Self::new(coeffs)
    }
}

// Tensor-product multiplication with numeric surface: result has degrees added in each direction
impl<T> Mul<BernsteinSurface<EFloat64>> for BernsteinSurface<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: BernsteinSurface<EFloat64>) -> Self::Output {
        let (du, dv) = (self.degree_u(), self.degree_v());
        let (eu, ev) = (rhs.degree_u(), rhs.degree_v());
        let mut coeffs = vec![vec![T::zero(); dv + ev + 1]; du + eu + 1];

        for k in 0..=du + eu {
            let i_min = k.saturating_sub(eu);
            let i_max = du.min(k);
            for l in 0..=dv + ev {
                let j_min = l.saturating_sub(ev);
                let j_max = dv.min(l);
                let mut acc = T::zero();
                for i in i_min..=i_max {
                    let a = k - i; // index in rhs u
                    for j in j_min..=j_max {
                        let b = l - j; // index in rhs v
                        let factor_u = (EFloat64::from(binomial_coefficient(du, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, a) as f64)
                            / EFloat64::from(binomial_coefficient(du + eu, k) as f64))
                        .unwrap();
                        let factor_v = (EFloat64::from(binomial_coefficient(dv, j) as f64)
                            * EFloat64::from(binomial_coefficient(ev, b) as f64)
                            / EFloat64::from(binomial_coefficient(dv + ev, l) as f64))
                        .unwrap();
                        let factor = factor_u * factor_v;
                        acc = acc
                            + self.coefficients[i][j].clone()
                                * (rhs.coefficients[a][b].clone() * factor);
                    }
                }
                coeffs[k][l] = acc;
            }
        }
        Self::new(coeffs)
    }
}

impl BernsteinSurface<Point> {
    pub fn dot(&self, rhs: &Self) -> BernsteinSurface<EFloat64> {
        let (du, dv) = (self.degree_u(), self.degree_v());
        let (eu, ev) = (rhs.degree_u(), rhs.degree_v());
        let mut coeffs = vec![vec![EFloat64::zero(); dv + ev + 1]; du + eu + 1];

        for k in 0..=du + eu {
            let i_min = k.saturating_sub(eu);
            let i_max = du.min(k);
            for l in 0..=dv + ev {
                let j_min = l.saturating_sub(ev);
                let j_max = dv.min(l);
                let mut acc = EFloat64::zero();
                for i in i_min..=i_max {
                    let a = k - i;
                    for j in j_min..=j_max {
                        let b = l - j;
                        let fu = (EFloat64::from(binomial_coefficient(du, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, a) as f64)
                            / EFloat64::from(binomial_coefficient(du + eu, k) as f64))
                        .unwrap();
                        let fv = (EFloat64::from(binomial_coefficient(dv, j) as f64)
                            * EFloat64::from(binomial_coefficient(ev, b) as f64)
                            / EFloat64::from(binomial_coefficient(dv + ev, l) as f64))
                        .unwrap();
                        let factor = fu * fv;
                        acc = acc + self.coefficients[i][j].dot(rhs.coefficients[a][b]) * factor;
                    }
                }
                coeffs[k][l] = acc;
            }
        }
        BernsteinSurface::new(coeffs)
    }

    pub fn cross(&self, rhs: &Self) -> Self {
        let (du, dv) = (self.degree_u(), self.degree_v());
        let (eu, ev) = (rhs.degree_u(), rhs.degree_v());
        let mut coeffs = vec![vec![Point::zero(); dv + ev + 1]; du + eu + 1];

        for k in 0..=du + eu {
            let i_min = k.saturating_sub(eu);
            let i_max = du.min(k);
            for l in 0..=dv + ev {
                let j_min = l.saturating_sub(ev);
                let j_max = dv.min(l);
                let mut acc = Point::zero();
                for i in i_min..=i_max {
                    let a = k - i;
                    for j in j_min..=j_max {
                        let b = l - j;
                        let fu = (EFloat64::from(binomial_coefficient(du, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, a) as f64)
                            / EFloat64::from(binomial_coefficient(du + eu, k) as f64))
                        .unwrap();
                        let fv = (EFloat64::from(binomial_coefficient(dv, j) as f64)
                            * EFloat64::from(binomial_coefficient(ev, b) as f64)
                            / EFloat64::from(binomial_coefficient(dv + ev, l) as f64))
                        .unwrap();
                        let factor = fu * fv;
                        acc = acc + self.coefficients[i][j].cross(rhs.coefficients[a][b]) * factor;
                    }
                }
                coeffs[k][l] = acc;
            }
        }
        Self::new(coeffs)
    }
}

// Exact division by a numeric Bernstein surface along tensor-product (lhs / rhs)
impl<T> Div<BernsteinSurface<EFloat64>> for BernsteinSurface<T>
where
    T: Zero
        + Clone
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<EFloat64, Output = T>
        + Div<EFloat64, Output = AlgebraResult<T>>
        + PartialEq,
{
    type Output = AlgebraResult<Self>;
    fn div(self, rhs: BernsteinSurface<EFloat64>) -> AlgebraResult<Self> {
        // Strategy: generalize 1D division to 2D via forward substitution over anti-diagonals (k,l)
        // Assume rhs(0,0) != 0 to have stable start.
        let du = self.degree_u();
        let dv = self.degree_v();
        let eu = rhs.degree_u();
        let ev = rhs.degree_v();
        if du < eu || dv < ev {
            return Err("Division degree mismatch: dividend degree < divisor degree".into());
        }

        let qu = du - eu;
        let qv = dv - ev;
        let mut q = vec![vec![T::zero(); qv + 1]; qu + 1];

        let r00 = rhs.coefficients[0][0];
        if r00 == 0.0 {
            return Err("Division by surface with zero at (u=0,v=0) is not supported".into());
        }

        for k in 0..=qu {
            for l in 0..=qv {
                // S = sum over (i,j)!=(k,l) with i<=k, j<=l of q[i][j] * r[k-i][l-j] * alpha
                let mut s_acc = T::zero();
                for i in 0..=k {
                    for j in 0..=l {
                        if i == k && j == l {
                            continue;
                        }
                        if k - i > eu || l - j > ev {
                            continue;
                        }
                        let alpha_u = (EFloat64::from(binomial_coefficient(qu, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, k - i) as f64)
                            / EFloat64::from(binomial_coefficient(qu + eu, k) as f64))
                        .unwrap();
                        let alpha_v = (EFloat64::from(binomial_coefficient(qv, j) as f64)
                            * EFloat64::from(binomial_coefficient(ev, l - j) as f64)
                            / EFloat64::from(binomial_coefficient(qv + ev, l) as f64))
                        .unwrap();
                        let alpha = alpha_u * alpha_v;
                        s_acc = s_acc + q[i][j].clone() * (rhs.coefficients[k - i][l - j] * alpha);
                    }
                }

                let alpha_k = (EFloat64::from(binomial_coefficient(qu, k) as f64)
                    / EFloat64::from(binomial_coefficient(qu + eu, k) as f64))
                .unwrap();
                let alpha_l = (EFloat64::from(binomial_coefficient(qv, l) as f64)
                    / EFloat64::from(binomial_coefficient(qv + ev, l) as f64))
                .unwrap();
                let denom = r00 * (alpha_k * alpha_l);

                let numer = self.coefficients[k][l].clone() - s_acc;
                q[k][l] = (numer / denom)?;
            }
        }

        let quotient = BernsteinSurface::new(q);
        let recomposed: BernsteinSurface<T> = quotient.clone() * rhs;
        if recomposed != self {
            return Err("Division has a remainder; not exactly divisible".into());
        }
        Ok(quotient)
    }
}

impl BernsteinSurface<EFloat64> {
    pub fn sqrt(&self) -> AlgebraResult<BernsteinSurface<EFloat64>> {
        // Require both degrees even and perfect square tensor-product wise
        let du = self.degree_u();
        let dv = self.degree_v();
        if du % 2 == 1 || dv % 2 == 1 {
            return Err(
                "Square root only defined for even degrees (perfect square) surfaces".into(),
            );
        }
        let mu = du / 2;
        let mv = dv / 2;

        let c = &self.coefficients; // c[i][j]
        if c[0][0] == 0.0 {
            // simplistic contradiction check similar to curve case
            if du >= 1 && c[1][0] != 0.0 {
                return Err("Surface not a perfect square (first non-zero index odd in u)".into());
            }
            if dv >= 1 && c[0][1] != 0.0 {
                return Err("Surface not a perfect square (first non-zero index odd in v)".into());
            }
        }

        let d00 = c[0][0].sqrt()?;
        let mut d = vec![vec![EFloat64::zero(); mv + 1]; mu + 1];
        d[0][0] = d00;

        // Solve for d[k,l] increasing in (k+l), excluding (0,0)
        for s in 0..=(mu + mv) {
            for k in 0..=s.min(mu) {
                let l = s - k;
                if l > mv {
                    continue;
                }
                if k == 0 && l == 0 {
                    continue;
                }

                // Accumulate S = sum over i=0..k, j=0..l, excluding (k,l):
                // factor(mu,i,k-i) * factor(mv,j,l-j) * d[i][j] * d[k-i][l-j]
                let mut s_acc = EFloat64::zero();
                for i in 0..=k {
                    for j in 0..=l {
                        if i == k && j == l {
                            continue;
                        }
                        let fu = (EFloat64::from(binomial_coefficient(mu, i) as f64)
                            * EFloat64::from(binomial_coefficient(mu, k - i) as f64)
                            / EFloat64::from(binomial_coefficient(2 * mu, k) as f64))
                        .unwrap();
                        let fv = (EFloat64::from(binomial_coefficient(mv, j) as f64)
                            * EFloat64::from(binomial_coefficient(mv, l - j) as f64)
                            / EFloat64::from(binomial_coefficient(2 * mv, l) as f64))
                        .unwrap();
                        s_acc = s_acc + (d[i][j] * d[k - i][l - j] * (fu * fv));
                    }
                }

                // Linear coefficient for d[k,l]: 2 * (C(mu,k)/C(2mu,k)) * (C(mv,l)/C(2mv,l)) * d00
                let au = (EFloat64::from(binomial_coefficient(mu, k) as f64)
                    / EFloat64::from(binomial_coefficient(2 * mu, k) as f64))
                .unwrap();
                let av = (EFloat64::from(binomial_coefficient(mv, l) as f64)
                    / EFloat64::from(binomial_coefficient(2 * mv, l) as f64))
                .unwrap();
                let denom = EFloat64::from(2.0) * au * av * d[0][0];

                let ckl = c[k][l];
                d[k][l] = ((ckl - s_acc) / denom)?;
            }
        }

        let q = BernsteinSurface::new(d);
        let back = q.clone() * q.clone();
        if back != *self {
            return Err("Verification of surface square root failed".into());
        }
        Ok(q)
    }
}

impl SurfaceLike for BernsteinSurface<Point> {
    fn u_span(&self) -> (EFloat64, EFloat64) {
        (EFloat64::zero(), EFloat64::one())
    }
    fn v_span(&self) -> (EFloat64, EFloat64) {
        (EFloat64::zero(), EFloat64::one())
    }
    fn eval(&self, u: EFloat64, v: EFloat64) -> Point {
        self.eval(u, v)
    }

    fn subdivide_u(&self, t: EFloat64) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized,
    {
        Ok(self.subdivide_u(t))
    }
    fn subdivide_v(&self, t: EFloat64) -> AlgebraResult<(Self, Self)>
    where
        Self: Sized,
    {
        Ok(self.subdivide_v(t))
    }

    fn get_convex_hull(&self) -> AlgebraResult<ConvexHull> {
        let points: Vec<Point> = self
            .coefficients
            .iter()
            .flat_map(|row| row.iter())
            .cloned()
            .collect();
        ConvexHull::try_new(points)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{color::Color10, line::Line, primitive_scene::PrimitiveScene};

    fn sample_surface() -> BernsteinSurface<Point> {
        // 4x4 grid (degree 3,3)
        let coeffs = vec![
            vec![
                Point::from_f64(0.0, 0.0, 0.0),
                Point::from_f64(0.0, 1.0, 0.0),
                Point::from_f64(0.0, 2.0, 0.0),
                Point::from_f64(0.0, 3.0, 0.0),
            ],
            vec![
                Point::from_f64(1.0, 0.0, 0.5),
                Point::from_f64(1.0, 1.0, 1.0),
                Point::from_f64(1.0, 2.0, 0.5),
                Point::from_f64(1.0, 3.0, 0.0),
            ],
            vec![
                Point::from_f64(2.0, 0.0, 0.0),
                Point::from_f64(2.0, 1.0, 0.5),
                Point::from_f64(2.0, 2.0, 1.0),
                Point::from_f64(2.0, 3.0, 0.5),
            ],
            vec![
                Point::from_f64(3.0, 0.0, 0.0),
                Point::from_f64(3.0, 1.0, 0.0),
                Point::from_f64(3.0, 2.0, 0.0),
                Point::from_f64(3.0, 3.0, 0.0),
            ],
        ];
        BernsteinSurface::new(coeffs)
    }

    #[test]
    fn test_eval_center() {
        let s = sample_surface();
        let p = s.eval(EFloat64::from(0.5), EFloat64::from(0.5));
        assert!(p != Point::zero());
    }

    #[test]
    fn test_elevate_and_eval_invariance() {
        let s = sample_surface();
        let s2 = s.elevate_degree(1, 2);
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                assert_eq!(s.eval(u.clone(), v.clone()), s2.eval(u, v));
            }
        }
    }

    #[test]
    fn test_subdivide_u_v() {
        let s = sample_surface();
        let t = EFloat64::from(0.5);
        let (lu, ru) = s.subdivide_u(t.clone());
        let (lv, rv) = s.subdivide_v(t.clone());
        // mapping: left-u covers [0,t] globally, so s(u*t,v) == left(u,v)
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                let su = s.eval(u.clone() * t.clone(), v.clone());
                let luu = lu.eval(u.clone(), v.clone());
                assert_eq!(su, luu);
            }
        }
        // mapping: right-u covers [t,1] globally, so s(t + u*(1-t), v) == right(u,v)
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                let one_minus_t = EFloat64::one() - t.clone();
                let mapped = t.clone() + (u.clone() * one_minus_t);
                let su = s.eval(mapped, v.clone());
                let ruu = ru.eval(u.clone(), v.clone());
                assert_eq!(su, ruu);
            }
        }

        // Subdivide in v: mapping analogous
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                let sv = s.eval(u.clone(), v.clone() * t.clone());
                let lvv = lv.eval(u.clone(), v.clone());
                assert_eq!(sv, lvv);
            }
        }
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                let one_minus_t = EFloat64::one() - t.clone();
                let mapped = t.clone() + (v.clone() * one_minus_t);
                let sv = s.eval(u.clone(), mapped);
                let rvv = rv.eval(u.clone(), v.clone());
                assert_eq!(sv, rvv);
            }
        }
    }

    #[test]
    fn test_derivatives_surfaces() -> AlgebraResult<()> {
        let s = sample_surface();
        let su = s.derivative_u();
        let sv = s.derivative_v();
        // visualize few tangents
        let mut scene = PrimitiveScene::new();
        scene.add_surface_like(&s, Color10::Blue, 30)?;
        for i in 0..=4 {
            for j in 0..=4 {
                let u = EFloat64::from(i as f64 / 4.0);
                let v = EFloat64::from(j as f64 / 4.0);
                let p = s.eval(u.clone(), v.clone());
                let tu = su.eval(u.clone(), v.clone()) * EFloat64::from(0.1);
                let tv = sv.eval(u.clone(), v.clone()) * EFloat64::from(0.1);
                scene.add_line(Line::try_new(p, p + tu).unwrap(), Color10::Red);
                scene.add_line(Line::try_new(p, p + tv).unwrap(), Color10::Green);
            }
        }
        scene.save_to_file("test_outputs/bernstein_surface_derivative.html")?;
        Ok(())
    }

    #[test]
    fn test_add_sub_scalar_mul() -> AlgebraResult<()> {
        let s = sample_surface();
        let s2 = sample_surface();
        let sum = s.clone() + s2.clone();
        let diff = s.clone() - s2.clone();
        let sm = s.clone() * EFloat64::from(2.5);
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                assert_eq!(
                    sum.eval(u.clone(), v.clone()),
                    s.eval(u.clone(), v.clone()) + s2.eval(u.clone(), v.clone())
                );
                assert_eq!(
                    diff.eval(u.clone(), v.clone()),
                    s.eval(u.clone(), v.clone()) - s2.eval(u.clone(), v.clone())
                );
                assert_eq!(
                    sm.eval(u.clone(), v.clone()),
                    s.eval(u.clone(), v.clone()) * EFloat64::from(2.5)
                );
            }
        }
        Ok(())
    }

    #[test]
    fn test_dot_cross_mul_div() -> AlgebraResult<()> {
        let s = sample_surface();
        let s2 = sample_surface();
        let dot = s.dot(&s2);
        let cross = s.cross(&s2);
        // product/division round trip for numeric case

        let q = BernsteinSurface::new(vec![vec![EFloat64::from(0.5), EFloat64::from(0.5)]]); // degree (0,1)
        let p = s.clone() * q.clone();
        let r = (p.clone() / q.clone())?;
        assert_eq!(r, s);

        // Evaluate dot product equality in a grid
        for iu in 0..=10 {
            for iv in 0..=10 {
                let u = EFloat64::from(iu as f64 / 10.0);
                let v = EFloat64::from(iv as f64 / 10.0);
                let s_dot_s2 = s
                    .eval(u.clone(), v.clone())
                    .dot(s2.eval(u.clone(), v.clone()));
                let dot_eval = dot.eval(u.clone(), v.clone());
                assert_eq!(s_dot_s2, dot_eval, "Dot product mismatch at ({}, {})", u, v);

                let s_cross_s2 = s
                    .eval(u.clone(), v.clone())
                    .cross(s2.eval(u.clone(), v.clone()));
                let cross_eval = cross.eval(u.clone(), v.clone());
                assert_eq!(
                    s_cross_s2, cross_eval,
                    "Cross product mismatch at ({}, {})",
                    u, v
                );
            }
        }

        Ok(())
    }

    // PH-like surface toy: derivative fields have polynomial speed; check sqrt on EFloat surface
    #[test]
    fn test_speed_sqrt_surface_like() -> AlgebraResult<()> {
        let s = sample_surface();
        let su = s.derivative_u();
        let sv = s.derivative_v();
        let e = su.dot(&su) + sv.dot(&sv); // squared speed of param mapping
        // Construct a square by multiplying a lower-degree tensor poly with itself
        let q = BernsteinSurface::new(vec![
            vec![EFloat64::from(0.5), EFloat64::from(1.0)],
            vec![EFloat64::from(1.0), EFloat64::from(2.0)],
        ]); // degree (1,1)
        let p = q.clone() * q.clone(); // degree (2,2)
        let srt = p.sqrt()?;
        let back = srt.clone() * srt.clone();
        assert_eq!(back, p);
        // just use e to touch the API
        let _ = e.eval(EFloat64::from(0.25), EFloat64::from(0.75));
        Ok(())
    }

    // Build a PH-like cubic surface from two complex hodographs for u and v.
    // Surface is constructed as separable sum: S(u,v) = C_u(u) + C_v(v) - P00,
    // so that S_u depends only on u and S_v only on v, each with PH-like speed.
    fn ph_cubic_surface_from_w(
        p00: (f64, f64),
        wu0: (f64, f64),
        wu1: (f64, f64),
        wv0: (f64, f64),
        wv1: (f64, f64),
    ) -> AlgebraResult<BernsteinSurface<Point>> {
        let base = Point::from_f64(p00.0, p00.1, 0.0);

        let sq = |a: f64, b: f64| (a * a - b * b, 2.0 * a * b);
        let mul = |a: f64, b: f64, c: f64, d: f64| (a * c - b * d, a * d + b * c);

        // U-direction PH cubic curve control points
        let (a0u, b0u) = wu0;
        let (a1u, b1u) = wu1;
        let d0u = sq(a0u, b0u);
        let d1u = mul(a0u, b0u, a1u, b1u);
        let d2u = sq(a1u, b1u);
        let n = EFloat64::from(3.0);
        let cu0 = base;
        let cu1 = cu0 + (Point::from_f64(d0u.0, d0u.1, 0.0) / n)?;
        let cu2 = cu1 + (Point::from_f64(d1u.0, d1u.1, 0.0) / n)?;
        let cu3 = cu2 + (Point::from_f64(d2u.0, d2u.1, 0.0) / n)?;

        // V-direction PH cubic curve control points
        let (a0v, b0v) = wv0;
        let (a1v, b1v) = wv1;
        let d0v = sq(a0v, b0v);
        let d1v = mul(a0v, b0v, a1v, b1v);
        let d2v = sq(a1v, b1v);
        let cv0 = base;
        let cv1 = cv0 + (Point::from_f64(0.0, d0v.1, d0v.0) / n)?;
        let cv2 = cv1 + (Point::from_f64(0.0, d1v.1, d1v.0) / n)?;
        let cv3 = cv2 + (Point::from_f64(0.0, d2v.1, d2v.0) / n)?;

        let cu = [cu0, cu1, cu2, cu3];
        let cv = [cv0, cv1, cv2, cv3];

        // Tensor product grid by separable sum
        let mut coeffs = vec![vec![Point::zero(); 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                coeffs[i][j] = cu[i] + cv[j] - base;
            }
        }

        Ok(BernsteinSurface::new(coeffs))
    }

    #[test]
    fn test_ph_like_u_derivative_sqrt_and_render() -> AlgebraResult<()> {
        // Build PH-like cubic surface by separable construction in u and v
        let s =
            ph_cubic_surface_from_w((0.0, 0.0), (3.0, 2.0), (2.0, -1.0), (0.5, 1.0), (1.0, 0.5))?;
        let su = s.derivative_u(); // degree (2,3)
        let speed_u_sq = su.dot(&su); // degree (4,6)
        let speed_u_poly = speed_u_sq.sqrt().ok();

        // Render: surface wireframe and normalized u-derivatives using sqrt factor
        let mut scene = PrimitiveScene::new();
        scene.add_surface_like_wireframe(&s, Color10::Blue, 24)?;
        for iu in 0..=24 {
            for iv in 0..=24 {
                let u = EFloat64::from(iu as f64 / 24.0);
                let v = EFloat64::from(iv as f64 / 24.0);
                let base = s.eval(u.clone(), v.clone());
                let du_vec = su.eval(u.clone(), v.clone());
                let denom = match &speed_u_poly {
                    Some(p) => p.eval(u.clone(), v.clone()),
                    None => du_vec.norm(),
                };
                let tangent = if denom != EFloat64::zero() {
                    (du_vec / denom)?
                } else {
                    Point::zero()
                };
                let tip = base + tangent * EFloat64::from(0.25);
                if tip != base {
                    scene.add_line(Line::try_new(base, tip)?, Color10::Red);
                }
            }
        }
        scene.save_to_file("test_outputs/bernstein_surface_ph_u_like.html")?;
        Ok(())
    }

    #[test]
    fn test_surface_normals_render_with_sqrt() -> AlgebraResult<()> {
        let s = sample_surface();
        let su = s.derivative_u();
        let sv = s.derivative_v();
        let mut scene = PrimitiveScene::new();
        scene.add_surface_like(&s, Color10::Blue, 22)?;
        // Normalize tangents via pointwise EFloat sqrt (not polynomial sqrt)
        let su_sq = su.dot(&su);
        let sv_sq = sv.dot(&sv);
        for iu in 0..=22 {
            for iv in 0..=22 {
                let u = EFloat64::from(iu as f64 / 22.0);
                let v = EFloat64::from(iv as f64 / 22.0);
                let p = s.eval(u.clone(), v.clone());
                let uvec = su.eval(u.clone(), v.clone());
                let vvec = sv.eval(u.clone(), v.clone());
                let ul = su_sq.eval(u.clone(), v.clone()).sqrt()?;
                let vl = sv_sq.eval(u.clone(), v.clone()).sqrt()?;
                let tu = if ul != EFloat64::zero() {
                    (uvec / ul)?
                } else {
                    Point::zero()
                };
                let tv = if vl != EFloat64::zero() {
                    (vvec / vl)?
                } else {
                    Point::zero()
                };
                let n = tu.cross(tv);
                let tip_u = p + tu * EFloat64::from(0.2);
                let tip_v = p + tv * EFloat64::from(0.2);
                let tip_n = p + n * EFloat64::from(0.2);
                if tip_u != p {
                    scene.add_line(Line::try_new(p, tip_u)?, Color10::Red);
                }
                if tip_v != p {
                    scene.add_line(Line::try_new(p, tip_v)?, Color10::Green);
                }
                if tip_n != p {
                    scene.add_line(Line::try_new(p, tip_n)?, Color10::Purple);
                }
            }
        }
        scene.save_to_file("test_outputs/bernstein_surface_speed_normals.html")?;
        Ok(())
    }
}
