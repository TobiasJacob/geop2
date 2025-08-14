use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

use crate::algebra_error::AlgebraResult;
use crate::bernstein::bernstein_surface::BernsteinSurface;
use crate::binomial_coefficient;
use crate::primitives::{efloat::EFloat64, point::Point};
use crate::zero::Zero;

use super::bernstein_curve::BernsteinCurve;

/// Tensor-product Bernstein volume with coefficients laid out as a 3D grid.
/// Dimensions are (degree_u+1) x (degree_v+1) x (degree_w+1).
#[derive(Debug, Clone, PartialEq)]
pub struct BernsteinVolume<T> {
    pub coefficients: Vec<Vec<Vec<T>>>, // [i_u][i_v][i_w]
}

impl<T> BernsteinVolume<T> {
    pub fn new(coefficients: Vec<Vec<Vec<T>>>) -> Self {
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

    pub fn degree_w(&self) -> usize {
        if self.coefficients.is_empty() || self.coefficients[0].is_empty() {
            0
        } else {
            self.coefficients[0][0].len().saturating_sub(1)
        }
    }
}

impl<T> Display for BernsteinVolume<T>
where
    T: Display + Zero + PartialEq,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        write!(f, "BernsteinVolume(\n")?;
        for i in 0..=du {
            write!(f, "[")?;
            for j in 0..=dv {
                write!(f, "[")?;
                for k in 0..=dw {
                    let c = &self.coefficients[i][j][k];
                    write!(f, "{}", c)?;
                    if k < dw {
                        write!(f, ", ")?;
                    }
                }
                write!(f, "]")?;
                if j < dv {
                    write!(f, "\n")?;
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

impl<T> BernsteinVolume<T>
where
    T: Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    /// Evaluate via tensor de Casteljau: reduce along u, then v, then w.
    pub fn eval(&self, u: EFloat64, v: EFloat64, w: EFloat64) -> T {
        if self.coefficients.is_empty()
            || self.coefficients[0].is_empty()
            || self.coefficients[0][0].is_empty()
        {
            panic!("Empty BernsteinVolume");
        }
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();

        // reduce along u for each (j,k) to a single value per j,k
        let mut uv: Vec<Vec<T>> = Vec::with_capacity(nv);
        for j in 0..nv {
            let mut row: Vec<T> = Vec::with_capacity(nw);
            for k in 0..nw {
                let mut beta: Vec<T> = (0..nu)
                    .map(|i| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nu {
                    for i in 0..(nu - r) {
                        beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                            + beta[i + 1].clone() * u.clone();
                    }
                }
                row.push(beta[0].clone());
            }
            uv.push(row);
        }

        // reduce along v for each k to a single value per k
        let mut gamma: Vec<T> = Vec::with_capacity(nw);
        for k in 0..nw {
            let mut beta: Vec<T> = (0..nv).map(|j| uv[j][k].clone()).collect();
            for r in 1..nv {
                for j in 0..(nv - r) {
                    beta[j] = beta[j].clone() * (EFloat64::one() - v.clone())
                        + beta[j + 1].clone() * v.clone();
                }
            }
            gamma.push(beta[0].clone());
        }

        // reduce along w to a single value
        let mut delta = gamma;
        for r in 1..nw {
            for k in 0..(nw - r) {
                delta[k] = delta[k].clone() * (EFloat64::one() - w.clone())
                    + delta[k + 1].clone() * w.clone();
            }
        }
        delta[0].clone()
    }

    /// Extract isoparametric curve in u for fixed (v,w)
    pub fn iso_u_at(&self, v: EFloat64, w: EFloat64) -> BernsteinCurve<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();

        // For each i along u, reduce the (v,w) plane to a single control point
        let mut cps: Vec<T> = Vec::with_capacity(nu);
        for i in 0..nu {
            // reduce along v for fixed i and each k
            let mut vw_reduced: Vec<T> = Vec::with_capacity(nw);
            for k in 0..nw {
                let mut beta: Vec<T> = (0..nv)
                    .map(|j| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nv {
                    for j in 0..(nv - r) {
                        beta[j] = beta[j].clone() * (EFloat64::one() - v.clone())
                            + beta[j + 1].clone() * v.clone();
                    }
                }
                vw_reduced.push(beta[0].clone());
            }
            // reduce along w on vw_reduced
            let mut delta = vw_reduced;
            for r in 1..nw {
                for k in 0..(nw - r) {
                    delta[k] = delta[k].clone() * (EFloat64::one() - w.clone())
                        + delta[k + 1].clone() * w.clone();
                }
            }
            cps.push(delta[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    /// Extract isoparametric curve in v for fixed (u,w)
    pub fn iso_v_at(&self, u: EFloat64, w: EFloat64) -> BernsteinCurve<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let mut cps: Vec<T> = Vec::with_capacity(nv);
        for j in 0..nv {
            // reduce along u for fixed j and each k
            let mut uw_reduced: Vec<T> = Vec::with_capacity(nw);
            for k in 0..nw {
                let mut beta: Vec<T> = (0..nu)
                    .map(|i| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nu {
                    for i in 0..(nu - r) {
                        beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                            + beta[i + 1].clone() * u.clone();
                    }
                }
                uw_reduced.push(beta[0].clone());
            }
            let mut delta = uw_reduced;
            for r in 1..nw {
                for k in 0..(nw - r) {
                    delta[k] = delta[k].clone() * (EFloat64::one() - w.clone())
                        + delta[k + 1].clone() * w.clone();
                }
            }
            cps.push(delta[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    /// Extract isoparametric curve in w for fixed (u,v)
    pub fn iso_w_at(&self, u: EFloat64, v: EFloat64) -> BernsteinCurve<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let mut cps: Vec<T> = Vec::with_capacity(nw);
        for k in 0..nw {
            // reduce along u for fixed (j,k)
            let mut col_u: Vec<T> = Vec::with_capacity(nv);
            for j in 0..nv {
                let mut beta: Vec<T> = (0..nu)
                    .map(|i| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nu {
                    for i in 0..(nu - r) {
                        beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                            + beta[i + 1].clone() * u.clone();
                    }
                }
                col_u.push(beta[0].clone());
            }
            // reduce along v on col_u
            let mut delta = col_u;
            for r in 1..nv {
                for j in 0..(nv - r) {
                    delta[j] = delta[j].clone() * (EFloat64::one() - v.clone())
                        + delta[j + 1].clone() * v.clone();
                }
            }
            cps.push(delta[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    /// Extract isoparametric surface at fixed w
    pub fn iso_uv_at(&self, w: EFloat64) -> BernsteinSurface<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let mut coeffs: Vec<Vec<T>> = vec![vec![self.coefficients[0][0][0].clone(); nv]; nu];
        for i in 0..nu {
            for j in 0..nv {
                let mut beta: Vec<T> = (0..nw)
                    .map(|k| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nw {
                    for k in 0..(nw - r) {
                        beta[k] = beta[k].clone() * (EFloat64::one() - w.clone())
                            + beta[k + 1].clone() * w.clone();
                    }
                }
                coeffs[i][j] = beta[0].clone();
            }
        }
        BernsteinSurface::new(coeffs)
    }

    /// Extract isoparametric surface at fixed v
    pub fn iso_uw_at(&self, v: EFloat64) -> BernsteinSurface<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let mut coeffs: Vec<Vec<T>> = vec![vec![self.coefficients[0][0][0].clone(); nw]; nu];
        for i in 0..nu {
            for k in 0..nw {
                let mut beta: Vec<T> = (0..nv)
                    .map(|j| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nv {
                    for j in 0..(nv - r) {
                        beta[j] = beta[j].clone() * (EFloat64::one() - v.clone())
                            + beta[j + 1].clone() * v.clone();
                    }
                }
                coeffs[i][k] = beta[0].clone();
            }
        }
        BernsteinSurface::new(coeffs)
    }

    /// Extract isoparametric surface at fixed u
    pub fn iso_vw_at(&self, u: EFloat64) -> BernsteinSurface<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let mut coeffs: Vec<Vec<T>> = vec![vec![self.coefficients[0][0][0].clone(); nw]; nv];
        for j in 0..nv {
            for k in 0..nw {
                let mut beta: Vec<T> = (0..nu)
                    .map(|i| self.coefficients[i][j][k].clone())
                    .collect();
                for r in 1..nu {
                    for i in 0..(nu - r) {
                        beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                            + beta[i + 1].clone() * u.clone();
                    }
                }
                coeffs[j][k] = beta[0].clone();
            }
        }
        BernsteinSurface::new(coeffs)
    }
}

impl<T> BernsteinVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    /// Degree elevation by r_u, r_v, r_w in the respective directions
    pub fn elevate_degree(&self, r_u: usize, r_v: usize, r_w: usize) -> Self {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();

        // Elevate in u
        let mut elev_u: Vec<Vec<Vec<T>>> =
            vec![vec![vec![T::zero(); dw + 1]; dv + 1]; du + r_u + 1];
        for i_new in 0..=du + r_u {
            for i_old in i_new.saturating_sub(r_u)..=du.min(i_new) {
                let alpha = (EFloat64::from(
                    (binomial_coefficient(r_u, i_new - i_old) * binomial_coefficient(du, i_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(du + r_u, i_new) as f64))
                .unwrap();
                for j in 0..=dv {
                    for k in 0..=dw {
                        elev_u[i_new][j][k] = elev_u[i_new][j][k].clone()
                            + self.coefficients[i_old][j][k].clone() * alpha.clone();
                    }
                }
            }
        }

        // Elevate in v
        let mut elev_v: Vec<Vec<Vec<T>>> =
            vec![vec![vec![T::zero(); dw + 1]; dv + r_v + 1]; du + r_u + 1];
        for j_new in 0..=dv + r_v {
            for j_old in j_new.saturating_sub(r_v)..=dv.min(j_new) {
                let beta = (EFloat64::from(
                    (binomial_coefficient(r_v, j_new - j_old) * binomial_coefficient(dv, j_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(dv + r_v, j_new) as f64))
                .unwrap();
                for i in 0..=du + r_u {
                    for k in 0..=dw {
                        elev_v[i][j_new][k] = elev_v[i][j_new][k].clone()
                            + elev_u[i][j_old][k].clone() * beta.clone();
                    }
                }
            }
        }

        // Elevate in w
        let mut elev_w: Vec<Vec<Vec<T>>> =
            vec![vec![vec![T::zero(); dw + r_w + 1]; dv + r_v + 1]; du + r_u + 1];
        for k_new in 0..=dw + r_w {
            for k_old in k_new.saturating_sub(r_w)..=dw.min(k_new) {
                let gamma = (EFloat64::from(
                    (binomial_coefficient(r_w, k_new - k_old) * binomial_coefficient(dw, k_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(dw + r_w, k_new) as f64))
                .unwrap();
                for i in 0..=du + r_u {
                    for j in 0..=dv + r_v {
                        elev_w[i][j][k_new] = elev_w[i][j][k_new].clone()
                            + elev_v[i][j][k_old].clone() * gamma.clone();
                    }
                }
            }
        }

        Self::new(elev_w)
    }

    /// Subdivide in u at parameter t using tensor de Casteljau
    pub fn subdivide_u(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();

        let mut left: Vec<Vec<Vec<T>>> = vec![vec![vec![T::zero(); dw]; dv]; du];
        let mut right: Vec<Vec<Vec<T>>> = vec![vec![vec![T::zero(); dw]; dv]; du];

        for j in 0..dv {
            for k in 0..dw {
                let mut beta: Vec<T> = (0..du)
                    .map(|i| self.coefficients[i][j][k].clone())
                    .collect();
                left[0][j][k] = beta[0].clone();
                right[du - 1][j][k] = beta[du - 1].clone();
                for r in 1..du {
                    for i in 0..(du - r) {
                        beta[i] = beta[i].clone() * (EFloat64::one() - t.clone())
                            + beta[i + 1].clone() * t.clone();
                    }
                    left[r][j][k] = beta[0].clone();
                    right[du - r - 1][j][k] = beta[du - r - 1].clone();
                }
            }
        }

        (Self::new(left), Self::new(right))
    }

    /// Subdivide in v at parameter t
    pub fn subdivide_v(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();

        let mut left: Vec<Vec<Vec<T>>> = vec![vec![vec![T::zero(); dw]; dv]; du];
        let mut right: Vec<Vec<Vec<T>>> = vec![vec![vec![T::zero(); dw]; dv]; du];

        for i in 0..du {
            for k in 0..dw {
                let mut beta: Vec<T> = (0..dv)
                    .map(|j| self.coefficients[i][j][k].clone())
                    .collect();
                left[i][0][k] = beta[0].clone();
                right[i][dv - 1][k] = beta[dv - 1].clone();
                for r in 1..dv {
                    for j in 0..(dv - r) {
                        beta[j] = beta[j].clone() * (EFloat64::one() - t.clone())
                            + beta[j + 1].clone() * t.clone();
                    }
                    left[i][r][k] = beta[0].clone();
                    right[i][dv - r - 1][k] = beta[dv - r - 1].clone();
                }
            }
        }

        (Self::new(left), Self::new(right))
    }

    /// Subdivide in w at parameter t
    pub fn subdivide_w(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();

        let mut left: Vec<Vec<Vec<T>>> = vec![vec![vec![T::zero(); dw]; dv]; du];
        let mut right: Vec<Vec<Vec<T>>> = vec![vec![vec![T::zero(); dw]; dv]; du];

        for i in 0..du {
            for j in 0..dv {
                let mut beta: Vec<T> = (0..dw)
                    .map(|k| self.coefficients[i][j][k].clone())
                    .collect();
                left[i][j][0] = beta[0].clone();
                right[i][j][dw - 1] = beta[dw - 1].clone();
                for r in 1..dw {
                    for k in 0..(dw - r) {
                        beta[k] = beta[k].clone() * (EFloat64::one() - t.clone())
                            + beta[k + 1].clone() * t.clone();
                    }
                    left[i][j][r] = beta[0].clone();
                    right[i][j][dw - r - 1] = beta[dw - r - 1].clone();
                }
            }
        }

        (Self::new(left), Self::new(right))
    }
}

impl<T> BernsteinVolume<T>
where
    T: Zero + Clone + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn derivative_u(&self) -> BernsteinVolume<T> {
        let du = self.degree_u();
        if du == 0 {
            return BernsteinVolume::new(vec![self.coefficients[0].clone()]);
        }
        let dv = self.degree_v();
        let dw = self.degree_w();
        let scale = EFloat64::from(du as f64);
        let mut coeffs = vec![vec![vec![T::zero(); dw + 1]; dv + 1]; du];
        for i in 0..du {
            for j in 0..=dv {
                for k in 0..=dw {
                    let diff =
                        self.coefficients[i + 1][j][k].clone() - self.coefficients[i][j][k].clone();
                    coeffs[i][j][k] = diff * scale.clone();
                }
            }
        }
        BernsteinVolume::new(coeffs)
    }

    pub fn derivative_v(&self) -> BernsteinVolume<T> {
        let dv = self.degree_v();
        if dv == 0 {
            // Collapse v-dimension to length 1 for each (i, k)
            let du = self.degree_u();
            let dw = self.degree_w();
            let mut base = vec![vec![vec![T::zero(); dw + 1]; 1]; du + 1];
            for i in 0..=du {
                for k in 0..=dw {
                    base[i][0][k] = self.coefficients[i][0][k].clone();
                }
            }
            return BernsteinVolume::new(base);
        }
        let du = self.degree_u();
        let dw = self.degree_w();
        let scale = EFloat64::from(dv as f64);
        let mut coeffs = vec![vec![vec![T::zero(); dw + 1]; dv]; du + 1];
        for i in 0..=du {
            for j in 0..dv {
                for k in 0..=dw {
                    let diff =
                        self.coefficients[i][j + 1][k].clone() - self.coefficients[i][j][k].clone();
                    coeffs[i][j][k] = diff * scale.clone();
                }
            }
        }
        BernsteinVolume::new(coeffs)
    }

    pub fn derivative_w(&self) -> BernsteinVolume<T> {
        let dw = self.degree_w();
        if dw == 0 {
            // Collapse w-dimension to length 1 for each (i, j)
            let du = self.degree_u();
            let dv = self.degree_v();
            let mut base = vec![vec![vec![T::zero(); 1]; dv + 1]; du + 1];
            for i in 0..=du {
                for j in 0..=dv {
                    base[i][j][0] = self.coefficients[i][j][0].clone();
                }
            }
            return BernsteinVolume::new(base);
        }
        let du = self.degree_u();
        let dv = self.degree_v();
        let scale = EFloat64::from(dw as f64);
        let mut coeffs = vec![vec![vec![T::zero(); dw]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..dw {
                    let diff =
                        self.coefficients[i][j][k + 1].clone() - self.coefficients[i][j][k].clone();
                    coeffs[i][j][k] = diff * scale.clone();
                }
            }
        }
        BernsteinVolume::new(coeffs)
    }
}

impl<T> BernsteinVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    fn equalize_degree(
        mut self,
        mut rhs: BernsteinVolume<T>,
    ) -> (BernsteinVolume<T>, BernsteinVolume<T>) {
        let (du1, dv1, dw1) = (self.degree_u(), self.degree_v(), self.degree_w());
        let (du2, dv2, dw2) = (rhs.degree_u(), rhs.degree_v(), rhs.degree_w());
        if du1 < du2 {
            self = self.elevate_degree(du2 - du1, 0, 0);
        }
        if dv1 < dv2 {
            self = self.elevate_degree(0, dv2 - dv1, 0);
        }
        if dw1 < dw2 {
            self = self.elevate_degree(0, 0, dw2 - dw1);
        }
        let (du1, dv1, dw1) = (self.degree_u(), self.degree_v(), self.degree_w());
        if du2 < du1 {
            rhs = rhs.elevate_degree(du1 - du2, 0, 0);
        }
        if dv2 < dv1 {
            rhs = rhs.elevate_degree(0, dv1 - dv2, 0);
        }
        if dw2 < dw1 {
            rhs = rhs.elevate_degree(0, 0, dw1 - dw2);
        }
        (self, rhs)
    }
}

impl<T> Add for BernsteinVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let du = lhs.degree_u();
        let dv = lhs.degree_v();
        let dw = lhs.degree_w();
        let mut coeffs = vec![vec![vec![T::zero(); dw + 1]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..=dw {
                    coeffs[i][j][k] =
                        lhs.coefficients[i][j][k].clone() + rhs.coefficients[i][j][k].clone();
                }
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Sub for BernsteinVolume<T>
where
    T: Zero + Clone + Sub<Output = T> + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let du = lhs.degree_u();
        let dv = lhs.degree_v();
        let dw = lhs.degree_w();
        let mut coeffs = vec![vec![vec![T::zero(); dw + 1]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..=dw {
                    coeffs[i][j][k] =
                        lhs.coefficients[i][j][k].clone() - rhs.coefficients[i][j][k].clone();
                }
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Mul<EFloat64> for BernsteinVolume<T>
where
    T: Clone + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: EFloat64) -> Self::Output {
        let mut coeffs = self.coefficients.clone();
        for plane in &mut coeffs {
            for row in plane {
                for c in row {
                    *c = c.clone() * rhs;
                }
            }
        }
        Self::new(coeffs)
    }
}

// Tensor-product multiplication with numeric volume: result degrees add in each direction
impl<T> Mul<BernsteinVolume<EFloat64>> for BernsteinVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: BernsteinVolume<EFloat64>) -> Self::Output {
        let (du, dv, dw) = (self.degree_u(), self.degree_v(), self.degree_w());
        let (eu, ev, ew) = (rhs.degree_u(), rhs.degree_v(), rhs.degree_w());
        let mut coeffs = vec![vec![vec![T::zero(); dw + ew + 1]; dv + ev + 1]; du + eu + 1];

        for p in 0..=du + eu {
            let i_min = p.saturating_sub(eu);
            let i_max = du.min(p);
            for q in 0..=dv + ev {
                let j_min = q.saturating_sub(ev);
                let j_max = dv.min(q);
                for r in 0..=dw + ew {
                    let k_min = r.saturating_sub(ew);
                    let k_max = dw.min(r);
                    let mut acc = T::zero();
                    for i in i_min..=i_max {
                        let a = p - i;
                        let fu = (EFloat64::from(binomial_coefficient(du, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, a) as f64)
                            / EFloat64::from(binomial_coefficient(du + eu, p) as f64))
                        .unwrap();
                        for j in j_min..=j_max {
                            let b = q - j;
                            let fv = (EFloat64::from(binomial_coefficient(dv, j) as f64)
                                * EFloat64::from(binomial_coefficient(ev, b) as f64)
                                / EFloat64::from(binomial_coefficient(dv + ev, q) as f64))
                            .unwrap();
                            for k in k_min..=k_max {
                                let c = r - k;
                                let fw = (EFloat64::from(binomial_coefficient(dw, k) as f64)
                                    * EFloat64::from(binomial_coefficient(ew, c) as f64)
                                    / EFloat64::from(binomial_coefficient(dw + ew, r) as f64))
                                .unwrap();
                                let factor = fu * fv * fw;
                                acc = acc
                                    + self.coefficients[i][j][k].clone()
                                        * (rhs.coefficients[a][b][c].clone() * factor);
                            }
                        }
                    }
                    coeffs[p][q][r] = acc;
                }
            }
        }
        Self::new(coeffs)
    }
}

impl BernsteinVolume<Point> {
    pub fn dot(&self, rhs: &Self) -> BernsteinVolume<EFloat64> {
        let (du, dv, dw) = (self.degree_u(), self.degree_v(), self.degree_w());
        let (eu, ev, ew) = (rhs.degree_u(), rhs.degree_v(), rhs.degree_w());
        let mut coeffs = vec![vec![vec![EFloat64::zero(); dw + ew + 1]; dv + ev + 1]; du + eu + 1];

        for p in 0..=du + eu {
            let i_min = p.saturating_sub(eu);
            let i_max = du.min(p);
            for q in 0..=dv + ev {
                let j_min = q.saturating_sub(ev);
                let j_max = dv.min(q);
                for r in 0..=dw + ew {
                    let k_min = r.saturating_sub(ew);
                    let k_max = dw.min(r);
                    let mut acc = EFloat64::zero();
                    for i in i_min..=i_max {
                        let a = p - i;
                        let fu = (EFloat64::from(binomial_coefficient(du, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, a) as f64)
                            / EFloat64::from(binomial_coefficient(du + eu, p) as f64))
                        .unwrap();
                        for j in j_min..=j_max {
                            let b = q - j;
                            let fv = (EFloat64::from(binomial_coefficient(dv, j) as f64)
                                * EFloat64::from(binomial_coefficient(ev, b) as f64)
                                / EFloat64::from(binomial_coefficient(dv + ev, q) as f64))
                            .unwrap();
                            for k in k_min..=k_max {
                                let c = r - k;
                                let fw = (EFloat64::from(binomial_coefficient(dw, k) as f64)
                                    * EFloat64::from(binomial_coefficient(ew, c) as f64)
                                    / EFloat64::from(binomial_coefficient(dw + ew, r) as f64))
                                .unwrap();
                                let factor = fu * fv * fw;
                                acc = acc
                                    + self.coefficients[i][j][k].dot(rhs.coefficients[a][b][c])
                                        * factor;
                            }
                        }
                    }
                    coeffs[p][q][r] = acc;
                }
            }
        }
        BernsteinVolume::new(coeffs)
    }

    pub fn cross(&self, rhs: &Self) -> Self {
        let (du, dv, dw) = (self.degree_u(), self.degree_v(), self.degree_w());
        let (eu, ev, ew) = (rhs.degree_u(), rhs.degree_v(), rhs.degree_w());
        let mut coeffs = vec![vec![vec![Point::zero(); dw + ew + 1]; dv + ev + 1]; du + eu + 1];

        for p in 0..=du + eu {
            let i_min = p.saturating_sub(eu);
            let i_max = du.min(p);
            for q in 0..=dv + ev {
                let j_min = q.saturating_sub(ev);
                let j_max = dv.min(q);
                for r in 0..=dw + ew {
                    let k_min = r.saturating_sub(ew);
                    let k_max = dw.min(r);
                    let mut acc = Point::zero();
                    for i in i_min..=i_max {
                        let a = p - i;
                        let fu = (EFloat64::from(binomial_coefficient(du, i) as f64)
                            * EFloat64::from(binomial_coefficient(eu, a) as f64)
                            / EFloat64::from(binomial_coefficient(du + eu, p) as f64))
                        .unwrap();
                        for j in j_min..=j_max {
                            let b = q - j;
                            let fv = (EFloat64::from(binomial_coefficient(dv, j) as f64)
                                * EFloat64::from(binomial_coefficient(ev, b) as f64)
                                / EFloat64::from(binomial_coefficient(dv + ev, q) as f64))
                            .unwrap();
                            for k in k_min..=k_max {
                                let c = r - k;
                                let fw = (EFloat64::from(binomial_coefficient(dw, k) as f64)
                                    * EFloat64::from(binomial_coefficient(ew, c) as f64)
                                    / EFloat64::from(binomial_coefficient(dw + ew, r) as f64))
                                .unwrap();
                                let factor = fu * fv * fw;
                                acc = acc
                                    + self.coefficients[i][j][k].cross(rhs.coefficients[a][b][c])
                                        * factor;
                            }
                        }
                    }
                    coeffs[p][q][r] = acc;
                }
            }
        }
        Self::new(coeffs)
    }
}

// Exact division by a numeric Bernstein volume along tensor-product (lhs / rhs)
impl<T> Div<BernsteinVolume<EFloat64>> for BernsteinVolume<T>
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
    fn div(self, rhs: BernsteinVolume<EFloat64>) -> AlgebraResult<Self> {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        let eu = rhs.degree_u();
        let ev = rhs.degree_v();
        let ew = rhs.degree_w();
        if du < eu || dv < ev || dw < ew {
            return Err("Division degree mismatch: dividend degree < divisor degree".into());
        }

        let qu = du - eu;
        let qv = dv - ev;
        let qw = dw - ew;
        let mut q = vec![vec![vec![T::zero(); qw + 1]; qv + 1]; qu + 1];

        let r000 = rhs.coefficients[0][0][0];
        if r000 == 0.0 {
            return Err("Division by volume with zero at (0,0,0) is not supported".into());
        }

        for k in 0..=qu {
            for l in 0..=qv {
                for m in 0..=qw {
                    // accumulate S over all (i<=k, j<=l, s<=m) except (k,l,m)
                    let mut s_acc = T::zero();
                    for i in 0..=k {
                        for j in 0..=l {
                            for s in 0..=m {
                                if i == k && j == l && s == m {
                                    continue;
                                }
                                if k - i > eu || l - j > ev || m - s > ew {
                                    continue;
                                }
                                let au = (EFloat64::from(binomial_coefficient(qu, i) as f64)
                                    * EFloat64::from(binomial_coefficient(eu, k - i) as f64)
                                    / EFloat64::from(binomial_coefficient(qu + eu, k) as f64))
                                .unwrap();
                                let av = (EFloat64::from(binomial_coefficient(qv, j) as f64)
                                    * EFloat64::from(binomial_coefficient(ev, l - j) as f64)
                                    / EFloat64::from(binomial_coefficient(qv + ev, l) as f64))
                                .unwrap();
                                let aw = (EFloat64::from(binomial_coefficient(qw, s) as f64)
                                    * EFloat64::from(binomial_coefficient(ew, m - s) as f64)
                                    / EFloat64::from(binomial_coefficient(qw + ew, m) as f64))
                                .unwrap();
                                let alpha = au * av * aw;
                                s_acc = s_acc
                                    + q[i][j][s].clone()
                                        * (rhs.coefficients[k - i][l - j][m - s] * alpha);
                            }
                        }
                    }

                    let ak = (EFloat64::from(binomial_coefficient(qu, k) as f64)
                        / EFloat64::from(binomial_coefficient(qu + eu, k) as f64))
                    .unwrap();
                    let al = (EFloat64::from(binomial_coefficient(qv, l) as f64)
                        / EFloat64::from(binomial_coefficient(qv + ev, l) as f64))
                    .unwrap();
                    let am = (EFloat64::from(binomial_coefficient(qw, m) as f64)
                        / EFloat64::from(binomial_coefficient(qw + ew, m) as f64))
                    .unwrap();
                    let denom = r000 * (ak * al * am);

                    let numer = self.coefficients[k][l][m].clone() - s_acc;
                    q[k][l][m] = (numer / denom)?;
                }
            }
        }

        let quotient = BernsteinVolume::new(q);
        let recomposed: BernsteinVolume<T> = quotient.clone() * rhs;
        if recomposed != self {
            return Err("Division has a remainder; not exactly divisible".into());
        }
        Ok(quotient)
    }
}

impl BernsteinVolume<EFloat64> {
    pub fn sqrt(&self) -> AlgebraResult<BernsteinVolume<EFloat64>> {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        if du % 2 == 1 || dv % 2 == 1 || dw % 2 == 1 {
            return Err(
                "Square root only defined for even degrees (perfect square) volumes".into(),
            );
        }
        let mu = du / 2;
        let mv = dv / 2;
        let mw = dw / 2;

        let c = &self.coefficients; // c[i][j][k]
        if c[0][0][0] == 0.0 {
            if du >= 1 && c[1][0][0] != 0.0 {
                return Err("Volume not a perfect square (first non-zero index odd in u)".into());
            }
            if dv >= 1 && c[0][1][0] != 0.0 {
                return Err("Volume not a perfect square (first non-zero index odd in v)".into());
            }
            if dw >= 1 && c[0][0][1] != 0.0 {
                return Err("Volume not a perfect square (first non-zero index odd in w)".into());
            }
        }

        let d000 = c[0][0][0].sqrt()?;
        let mut d = vec![vec![vec![EFloat64::zero(); mw + 1]; mv + 1]; mu + 1];
        d[0][0][0] = d000;

        // Solve for d[a][b][c] in increasing a+b+c
        for s in 0..=(mu + mv + mw) {
            for a in 0..=s.min(mu) {
                for b in 0..=(s - a).min(mv) {
                    let cidx = s - a - b;
                    if cidx > mw {
                        continue;
                    }
                    if a == 0 && b == 0 && cidx == 0 {
                        continue;
                    }

                    let mut s_acc = EFloat64::zero();
                    for i in 0..=a {
                        for j in 0..=b {
                            for k in 0..=cidx {
                                if i == a && j == b && k == cidx {
                                    continue;
                                }
                                let fu = (EFloat64::from(binomial_coefficient(mu, i) as f64)
                                    * EFloat64::from(binomial_coefficient(mu, a - i) as f64)
                                    / EFloat64::from(binomial_coefficient(2 * mu, a) as f64))
                                .unwrap();
                                let fv = (EFloat64::from(binomial_coefficient(mv, j) as f64)
                                    * EFloat64::from(binomial_coefficient(mv, b - j) as f64)
                                    / EFloat64::from(binomial_coefficient(2 * mv, b) as f64))
                                .unwrap();
                                let fw = (EFloat64::from(binomial_coefficient(mw, k) as f64)
                                    * EFloat64::from(binomial_coefficient(mw, cidx - k) as f64)
                                    / EFloat64::from(binomial_coefficient(2 * mw, cidx) as f64))
                                .unwrap();
                                s_acc = s_acc
                                    + (d[i][j][k] * d[a - i][b - j][cidx - k] * (fu * fv * fw));
                            }
                        }
                    }

                    let au = (EFloat64::from(binomial_coefficient(mu, a) as f64)
                        / EFloat64::from(binomial_coefficient(2 * mu, a) as f64))
                    .unwrap();
                    let av = (EFloat64::from(binomial_coefficient(mv, b) as f64)
                        / EFloat64::from(binomial_coefficient(2 * mv, b) as f64))
                    .unwrap();
                    let aw = (EFloat64::from(binomial_coefficient(mw, cidx) as f64)
                        / EFloat64::from(binomial_coefficient(2 * mw, cidx) as f64))
                    .unwrap();
                    let denom = EFloat64::from(2.0) * au * av * aw * d[0][0][0];

                    let cabc = c[a][b][cidx];
                    d[a][b][cidx] = ((cabc - s_acc) / denom)?;
                }
            }
        }

        let q = BernsteinVolume::new(d);
        let back = q.clone() * q.clone();
        if back != *self {
            return Err("Verification of volume square root failed".into());
        }
        Ok(q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_volume() -> BernsteinVolume<Point> {
        // 3x3x3 grid (degree 2,2,2)
        let mut coeffs = vec![vec![vec![Point::zero(); 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    coeffs[i][j][k] = Point::from_f64(i as f64, j as f64, k as f64)
                        + Point::from_f64(0.0, (i + j) as f64 * 0.1, (j + k) as f64 * 0.05);
                }
            }
        }

        coeffs[1][1][1] = Point::from_f64(1.0, 2.0, 3.0);
        coeffs[1][1][2] = Point::from_f64(4.0, 1.0, 2.0);

        BernsteinVolume::new(coeffs)
    }

    #[test]
    fn test_eval_center() {
        let v = sample_volume();
        let p = v.eval(
            EFloat64::from(0.5),
            EFloat64::from(0.5),
            EFloat64::from(0.5),
        );
        assert!(p != Point::zero());
    }

    #[test]
    fn test_elevate_and_eval_invariance() {
        let v = sample_volume();
        let v2 = v.elevate_degree(1, 1, 1);
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    assert_eq!(v.eval(u.clone(), vv.clone(), w.clone()), v2.eval(u, vv, w));
                }
            }
        }
    }

    #[test]
    fn test_subdivide_u_v_w() {
        let v = sample_volume();
        let t = EFloat64::from(0.5);
        let (lu, ru) = v.subdivide_u(t.clone());
        let (lv, rv) = v.subdivide_v(t.clone());
        let (lw, rw) = v.subdivide_w(t.clone());

        // u-left mapping: global maps (u*t, v, w)
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    let mapped = v.eval(u.clone() * t.clone(), vv.clone(), w.clone());
                    let sub = lu.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(mapped, sub);
                }
            }
        }
        // u-right mapping: t + u*(1-t)
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    let one_minus_t = EFloat64::one() - t.clone();
                    let mapped_u = t.clone() + (u.clone() * one_minus_t.clone());
                    let mapped = v.eval(mapped_u, vv.clone(), w.clone());
                    let sub = ru.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(mapped, sub);
                }
            }
        }

        // v-left mapping
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    let mapped = v.eval(u.clone(), vv.clone() * t.clone(), w.clone());
                    let sub = lv.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(mapped, sub);
                }
            }
        }
        // v-right mapping
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    let one_minus_t = EFloat64::one() - t.clone();
                    let mapped_v = t.clone() + (vv.clone() * one_minus_t.clone());
                    let mapped = v.eval(u.clone(), mapped_v, w.clone());
                    let sub = rv.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(mapped, sub);
                }
            }
        }

        // w-left mapping
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    let mapped = v.eval(u.clone(), vv.clone(), w.clone() * t.clone());
                    let sub = lw.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(mapped, sub);
                }
            }
        }
        // w-right mapping
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    let one_minus_t = EFloat64::one() - t.clone();
                    let mapped_w = t.clone() + (w.clone() * one_minus_t.clone());
                    let mapped = v.eval(u.clone(), vv.clone(), mapped_w);
                    let sub = rw.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(mapped, sub);
                }
            }
        }
    }

    #[test]
    fn test_derivatives() {
        let v = sample_volume();
        let vu = v.derivative_u();
        let vv = v.derivative_v();
        let vw = v.derivative_w();
        // Just evaluate at a few points to touch API
        for iu in 0..=2 {
            for iv in 0..=2 {
                for iw in 0..=2 {
                    let u = EFloat64::from(iu as f64 / 2.0);
                    let vpar = EFloat64::from(iv as f64 / 2.0);
                    let w = EFloat64::from(iw as f64 / 2.0);
                    let _ = vu.eval(u.clone(), vpar.clone(), w.clone());
                    let _ = vv.eval(u.clone(), vpar.clone(), w.clone());
                    let _ = vw.eval(u.clone(), vpar.clone(), w.clone());
                }
            }
        }
    }

    #[test]
    fn test_add_sub_scalar_mul() {
        let v1 = sample_volume();
        let v2 = sample_volume();
        let sum = v1.clone() + v2.clone();
        let diff = v1.clone() - v2.clone();
        let sm = v1.clone() * EFloat64::from(2.0);
        for iu in 0..=4 {
            for iv in 0..=4 {
                for iw in 0..=4 {
                    let u = EFloat64::from(iu as f64 / 4.0);
                    let vv = EFloat64::from(iv as f64 / 4.0);
                    let w = EFloat64::from(iw as f64 / 4.0);
                    assert_eq!(
                        sum.eval(u.clone(), vv.clone(), w.clone()),
                        v1.eval(u.clone(), vv.clone(), w.clone())
                            + v2.eval(u.clone(), vv.clone(), w.clone())
                    );
                    assert_eq!(
                        diff.eval(u.clone(), vv.clone(), w.clone()),
                        v1.eval(u.clone(), vv.clone(), w.clone())
                            - v2.eval(u.clone(), vv.clone(), w.clone())
                    );
                    assert_eq!(
                        sm.eval(u.clone(), vv.clone(), w.clone()),
                        v1.eval(u.clone(), vv.clone(), w.clone()) * EFloat64::from(2.0)
                    );
                }
            }
        }
    }

    #[test]
    fn test_dot_cross_mul_div_sqrt() -> AlgebraResult<()> {
        let v = sample_volume();
        let v2 = sample_volume();
        let dot = v.dot(&v2);
        let cross = v.cross(&v2);

        // product/division round trip for numeric case
        let q = BernsteinVolume::new(vec![vec![vec![EFloat64::from(0.5), EFloat64::from(0.5)]]]); // degree (0,0,1)
        let p = v.clone() * q.clone();
        let r = (p.clone() / q.clone())?;
        assert_eq!(r, v);

        // Evaluate dot/cross equality on a small grid
        for iu in 0..=3 {
            for iv in 0..=3 {
                for iw in 0..=3 {
                    let u = EFloat64::from(iu as f64 / 3.0);
                    let vv = EFloat64::from(iv as f64 / 3.0);
                    let w = EFloat64::from(iw as f64 / 3.0);
                    let lhs_dot = v.eval(u.clone(), vv.clone(), w.clone()).dot(v2.eval(
                        u.clone(),
                        vv.clone(),
                        w.clone(),
                    ));
                    let rhs_dot = dot.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(lhs_dot, rhs_dot);

                    let lhs_cross = v.eval(u.clone(), vv.clone(), w.clone()).cross(v2.eval(
                        u.clone(),
                        vv.clone(),
                        w.clone(),
                    ));
                    let rhs_cross = cross.eval(u.clone(), vv.clone(), w.clone());
                    assert_eq!(lhs_cross, rhs_cross);
                }
            }
        }

        // sqrt: build a perfect square from a (1,1,1) volume
        let _qq = BernsteinVolume::new(vec![
            vec![vec![EFloat64::from(0.5), EFloat64::from(1.0)]],
            vec![vec![EFloat64::from(1.0), EFloat64::from(2.0)]],
        ]); // degree (1,0,1)
        // Make it full cubic by adding a middle v dimension as (1)
        let qq = BernsteinVolume::new(vec![
            vec![vec![EFloat64::from(0.5), EFloat64::from(1.0)]],
            vec![vec![EFloat64::from(1.0), EFloat64::from(2.0)]],
        ]);
        let pp = qq.clone() * qq.clone();
        let srt = pp.sqrt()?;
        let back = srt.clone() * srt.clone();
        assert_eq!(back, pp);

        Ok(())
    }
}
