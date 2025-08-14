use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

use crate::algebra_error::AlgebraResult;
use crate::bernstein::bernstein_surface::BernsteinSurface;
use crate::binomial_coefficient;
use crate::primitives::convex_hull::ConvexHull;
use crate::primitives::{efloat::EFloat64, point::Point};
use crate::zero::Zero;

use super::bernstein_curve::BernsteinCurve;
use super::bernstein_volume::BernsteinVolume;

/// 4D Tensor-product Bernstein hypervolume with coefficients laid out as a 4D grid.
/// Dimensions are (degree_u+1) x (degree_v+1) x (degree_w+1) x (degree_x+1).
#[derive(Debug, Clone, PartialEq)]
pub struct BernsteinHyperVolume<T> {
    pub coefficients: Vec<Vec<Vec<Vec<T>>>>, // [i][j][k][l]
}

impl<T> BernsteinHyperVolume<T> {
    pub fn new(coefficients: Vec<Vec<Vec<Vec<T>>>>) -> Self {
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
    pub fn degree_x(&self) -> usize {
        if self.coefficients.is_empty()
            || self.coefficients[0].is_empty()
            || self.coefficients[0][0].is_empty()
        {
            0
        } else {
            self.coefficients[0][0][0].len().saturating_sub(1)
        }
    }
}

impl BernsteinHyperVolume<Point> {
    pub fn get_convex_hull(&self) -> AlgebraResult<ConvexHull> {
        let flattened = self
            .coefficients
            .iter()
            .flatten()
            .flatten()
            .flatten()
            .map(|x| *x)
            .collect::<Vec<_>>();
        ConvexHull::try_new(flattened)
    }
}

impl<T> Display for BernsteinHyperVolume<T>
where
    T: Display + Zero + PartialEq,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        let dx = self.degree_x();
        write!(f, "BernsteinHyperVolume(\n")?;
        for i in 0..=du {
            write!(f, "[")?;
            for j in 0..=dv {
                write!(f, "[")?;
                for k in 0..=dw {
                    write!(f, "[")?;
                    for l in 0..=dx {
                        write!(f, "{}", self.coefficients[i][j][k][l])?;
                        if l < dx {
                            write!(f, ", ")?;
                        }
                    }
                    write!(f, "]")?;
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

impl<T> BernsteinHyperVolume<T>
where
    T: Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn eval(&self, u: EFloat64, v: EFloat64, w: EFloat64, x: EFloat64) -> T {
        if self.coefficients.is_empty()
            || self.coefficients[0].is_empty()
            || self.coefficients[0][0].is_empty()
            || self.coefficients[0][0][0].is_empty()
        {
            panic!("Empty BernsteinHyperVolume");
        }
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();

        // Reduce along u -> 3D grid in (v,w,x)
        let mut vwx: Vec<Vec<Vec<T>>> = Vec::with_capacity(nv);
        for j in 0..nv {
            let mut wx: Vec<Vec<T>> = Vec::with_capacity(nw);
            for k in 0..nw {
                // Do u reduction per x index
                let mut row_x: Vec<T> = Vec::with_capacity(nx);
                for l in 0..nx {
                    let mut beta: Vec<T> = (0..nu)
                        .map(|i| self.coefficients[i][j][k][l].clone())
                        .collect();
                    for r in 1..nu {
                        for i in 0..(nu - r) {
                            beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                                + beta[i + 1].clone() * u.clone();
                        }
                    }
                    row_x.push(beta[0].clone());
                }
                wx.push(row_x);
            }
            vwx.push(wx);
        }

        // Reduce along v -> 2D grid in (w,x)
        let mut wx: Vec<Vec<T>> = Vec::with_capacity(nw);
        for k in 0..nw {
            let mut row_x: Vec<T> = Vec::with_capacity(nx);
            for l in 0..nx {
                let mut beta: Vec<T> = (0..nv).map(|j| vwx[j][k][l].clone()).collect();
                for r in 1..nv {
                    for j in 0..(nv - r) {
                        beta[j] = beta[j].clone() * (EFloat64::one() - v.clone())
                            + beta[j + 1].clone() * v.clone();
                    }
                }
                row_x.push(beta[0].clone());
            }
            wx.push(row_x);
        }

        // Reduce along w -> 1D in x
        let mut xline: Vec<T> = Vec::with_capacity(nx);
        for l in 0..nx {
            let mut beta: Vec<T> = (0..nw).map(|k| wx[k][l].clone()).collect();
            for r in 1..nw {
                for k in 0..(nw - r) {
                    beta[k] = beta[k].clone() * (EFloat64::one() - w.clone())
                        + beta[k + 1].clone() * w.clone();
                }
            }
            xline.push(beta[0].clone());
        }

        // Reduce along x -> scalar
        let mut delta = xline;
        for r in 1..nx {
            for l in 0..(nx - r) {
                delta[l] = delta[l].clone() * (EFloat64::one() - x.clone())
                    + delta[l + 1].clone() * x.clone();
            }
        }
        delta[0].clone()
    }

    // 1D isoparametric curves
    pub fn iso_u_at(&self, v: EFloat64, w: EFloat64, x: EFloat64) -> BernsteinCurve<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();
        let mut cps: Vec<T> = Vec::with_capacity(nu);
        for i in 0..nu {
            // Build a vector along w by first reducing along x
            let mut wline: Vec<T> = Vec::with_capacity(nw);
            for k in 0..nw {
                let mut xbeta: Vec<T> = (0..nx)
                    .map(|l| {
                        // v-reduced above into xline; recompute here directly
                        // compute v-reduced value at given (i,k,l)
                        let mut vbeta: Vec<T> = (0..nv)
                            .map(|j| self.coefficients[i][j][k][l].clone())
                            .collect();
                        for r in 1..nv {
                            for j in 0..(nv - r) {
                                vbeta[j] = vbeta[j].clone() * (EFloat64::one() - v.clone())
                                    + vbeta[j + 1].clone() * v.clone();
                            }
                        }
                        vbeta[0].clone()
                    })
                    .collect();
                // reduce along x
                for r in 1..nx {
                    for l in 0..(nx - r) {
                        xbeta[l] = xbeta[l].clone() * (EFloat64::one() - x.clone())
                            + xbeta[l + 1].clone() * x.clone();
                    }
                }
                wline.push(xbeta[0].clone());
            }
            // reduce along w to a single value
            for r in 1..nw {
                for k in 0..(nw - r) {
                    wline[k] = wline[k].clone() * (EFloat64::one() - w.clone())
                        + wline[k + 1].clone() * w.clone();
                }
            }
            cps.push(wline[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    pub fn iso_v_at(&self, u: EFloat64, w: EFloat64, x: EFloat64) -> BernsteinCurve<T> {
        // symmetric to iso_u_at
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();
        let mut cps: Vec<T> = Vec::with_capacity(nv);
        for j in 0..nv {
            let mut wline: Vec<T> = Vec::with_capacity(nw);
            for k in 0..nw {
                // reduce along u for fixed (j,k,*) to xline
                let mut xbeta: Vec<T> = Vec::with_capacity(nx);
                for l in 0..nx {
                    let mut beta: Vec<T> = (0..nu)
                        .map(|i| self.coefficients[i][j][k][l].clone())
                        .collect();
                    for r in 1..nu {
                        for i in 0..(nu - r) {
                            beta[i] = beta[i].clone() * (EFloat64::one() - u.clone())
                                + beta[i + 1].clone() * u.clone();
                        }
                    }
                    xbeta.push(beta[0].clone());
                }
                // reduce along x
                for r in 1..nx {
                    for l in 0..(nx - r) {
                        xbeta[l] = xbeta[l].clone() * (EFloat64::one() - x.clone())
                            + xbeta[l + 1].clone() * x.clone();
                    }
                }
                wline.push(xbeta[0].clone());
            }
            // reduce along w
            for r in 1..nw {
                for k in 0..(nw - r) {
                    wline[k] = wline[k].clone() * (EFloat64::one() - w.clone())
                        + wline[k + 1].clone() * w.clone();
                }
            }
            cps.push(wline[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    pub fn iso_w_at(&self, u: EFloat64, v: EFloat64, x: EFloat64) -> BernsteinCurve<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();
        let mut cps: Vec<T> = Vec::with_capacity(nw);
        for k in 0..nw {
            // reduce along u then v to xline
            let mut xbeta: Vec<T> = Vec::with_capacity(nx);
            for l in 0..nx {
                // reduce along u for all j, then along v
                let mut vbeta: Vec<T> = Vec::with_capacity(nv);
                for j in 0..nv {
                    let mut ubeta: Vec<T> = (0..nu)
                        .map(|i| self.coefficients[i][j][k][l].clone())
                        .collect();
                    for r in 1..nu {
                        for i in 0..(nu - r) {
                            ubeta[i] = ubeta[i].clone() * (EFloat64::one() - u.clone())
                                + ubeta[i + 1].clone() * u.clone();
                        }
                    }
                    vbeta.push(ubeta[0].clone());
                }
                for r in 1..nv {
                    for j in 0..(nv - r) {
                        vbeta[j] = vbeta[j].clone() * (EFloat64::one() - v.clone())
                            + vbeta[j + 1].clone() * v.clone();
                    }
                }
                xbeta.push(vbeta[0].clone());
            }
            // reduce along x
            for r in 1..nx {
                for l in 0..(nx - r) {
                    xbeta[l] = xbeta[l].clone() * (EFloat64::one() - x.clone())
                        + xbeta[l + 1].clone() * x.clone();
                }
            }
            cps.push(xbeta[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    pub fn iso_x_at(&self, u: EFloat64, v: EFloat64, w: EFloat64) -> BernsteinCurve<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();
        let mut cps: Vec<T> = Vec::with_capacity(nx);
        for l in 0..nx {
            // reduce along u, then v, then w
            let mut wbeta: Vec<T> = Vec::with_capacity(nw);
            for k in 0..nw {
                let mut vbeta: Vec<T> = Vec::with_capacity(nv);
                for j in 0..nv {
                    let mut ubeta: Vec<T> = (0..nu)
                        .map(|i| self.coefficients[i][j][k][l].clone())
                        .collect();
                    for r in 1..nu {
                        for i in 0..(nu - r) {
                            ubeta[i] = ubeta[i].clone() * (EFloat64::one() - u.clone())
                                + ubeta[i + 1].clone() * u.clone();
                        }
                    }
                    vbeta.push(ubeta[0].clone());
                }
                for r in 1..nv {
                    for j in 0..(nv - r) {
                        vbeta[j] = vbeta[j].clone() * (EFloat64::one() - v.clone())
                            + vbeta[j + 1].clone() * v.clone();
                    }
                }
                wbeta.push(vbeta[0].clone());
            }
            for r in 1..nw {
                for k in 0..(nw - r) {
                    wbeta[k] = wbeta[k].clone() * (EFloat64::one() - w.clone())
                        + wbeta[k + 1].clone() * w.clone();
                }
            }
            cps.push(wbeta[0].clone());
        }
        BernsteinCurve::new(cps)
    }

    // 2D isoparametric surfaces for fixed two parameters (examples)
    pub fn iso_uv_at(&self, w: EFloat64, x: EFloat64) -> BernsteinSurface<T> {
        // reduce along w then x, leaving (u,v)
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();
        let mut coeffs = vec![vec![self.coefficients[0][0][0][0].clone(); nv]; nu];
        for i in 0..nu {
            for j in 0..nv {
                // reduce along w to xline
                let mut xbeta: Vec<T> = Vec::with_capacity(nx);
                for l in 0..nx {
                    let mut wbeta: Vec<T> = (0..nw)
                        .map(|k| self.coefficients[i][j][k][l].clone())
                        .collect();
                    for r in 1..nw {
                        for k in 0..(nw - r) {
                            wbeta[k] = wbeta[k].clone() * (EFloat64::one() - w.clone())
                                + wbeta[k + 1].clone() * w.clone();
                        }
                    }
                    xbeta.push(wbeta[0].clone());
                }
                // reduce along x
                for r in 1..nx {
                    for l in 0..(nx - r) {
                        xbeta[l] = xbeta[l].clone() * (EFloat64::one() - x.clone())
                            + xbeta[l + 1].clone() * x.clone();
                    }
                }
                coeffs[i][j] = xbeta[0].clone();
            }
        }
        BernsteinSurface::new(coeffs)
    }

    // 3D isoparametric hyper-slice: fix x to get a 3D BernsteinVolume
    pub fn iso_uvw_at(&self, x: EFloat64) -> BernsteinVolume<T> {
        let nu = self.coefficients.len();
        let nv = self.coefficients[0].len();
        let nw = self.coefficients[0][0].len();
        let nx = self.coefficients[0][0][0].len();
        let mut coeffs = vec![vec![vec![self.coefficients[0][0][0][0].clone(); nw]; nv]; nu];
        for i in 0..nu {
            for j in 0..nv {
                for k in 0..nw {
                    // reduce along x
                    let mut beta: Vec<T> = (0..nx)
                        .map(|l| self.coefficients[i][j][k][l].clone())
                        .collect();
                    for r in 1..nx {
                        for l in 0..(nx - r) {
                            beta[l] = beta[l].clone() * (EFloat64::one() - x.clone())
                                + beta[l + 1].clone() * x.clone();
                        }
                    }
                    coeffs[i][j][k] = beta[0].clone();
                }
            }
        }
        BernsteinVolume::new(coeffs)
    }
}

impl<T> BernsteinHyperVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn elevate_degree(&self, r_u: usize, r_v: usize, r_w: usize, r_x: usize) -> Self {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        let dx = self.degree_x();

        // elevate u
        let mut eu = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; dv + 1]; du + r_u + 1];
        for i_new in 0..=du + r_u {
            for i_old in i_new.saturating_sub(r_u)..=du.min(i_new) {
                let alpha = (EFloat64::from(
                    (binomial_coefficient(r_u, i_new - i_old) * binomial_coefficient(du, i_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(du + r_u, i_new) as f64))
                .unwrap();
                for j in 0..=dv {
                    for k in 0..=dw {
                        for l in 0..=dx {
                            eu[i_new][j][k][l] = eu[i_new][j][k][l].clone()
                                + self.coefficients[i_old][j][k][l].clone() * alpha.clone();
                        }
                    }
                }
            }
        }

        // elevate v
        let mut ev = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; dv + r_v + 1]; du + r_u + 1];
        for j_new in 0..=dv + r_v {
            for j_old in j_new.saturating_sub(r_v)..=dv.min(j_new) {
                let beta = (EFloat64::from(
                    (binomial_coefficient(r_v, j_new - j_old) * binomial_coefficient(dv, j_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(dv + r_v, j_new) as f64))
                .unwrap();
                for i in 0..=du + r_u {
                    for k in 0..=dw {
                        for l in 0..=dx {
                            ev[i][j_new][k][l] = ev[i][j_new][k][l].clone()
                                + eu[i][j_old][k][l].clone() * beta.clone();
                        }
                    }
                }
            }
        }

        // elevate w
        let mut ew =
            vec![vec![vec![vec![T::zero(); dx + 1]; dw + r_w + 1]; dv + r_v + 1]; du + r_u + 1];
        for k_new in 0..=dw + r_w {
            for k_old in k_new.saturating_sub(r_w)..=dw.min(k_new) {
                let gamma = (EFloat64::from(
                    (binomial_coefficient(r_w, k_new - k_old) * binomial_coefficient(dw, k_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(dw + r_w, k_new) as f64))
                .unwrap();
                for i in 0..=du + r_u {
                    for j in 0..=dv + r_v {
                        for l in 0..=dx {
                            ew[i][j][k_new][l] = ew[i][j][k_new][l].clone()
                                + ev[i][j][k_old][l].clone() * gamma.clone();
                        }
                    }
                }
            }
        }

        // elevate x
        let mut ex = vec![
            vec![vec![vec![T::zero(); dx + r_x + 1]; dw + r_w + 1]; dv + r_v + 1];
            du + r_u + 1
        ];
        for l_new in 0..=dx + r_x {
            for l_old in l_new.saturating_sub(r_x)..=dx.min(l_new) {
                let delta = (EFloat64::from(
                    (binomial_coefficient(r_x, l_new - l_old) * binomial_coefficient(dx, l_old))
                        as f64,
                ) / EFloat64::from(binomial_coefficient(dx + r_x, l_new) as f64))
                .unwrap();
                for i in 0..=du + r_u {
                    for j in 0..=dv + r_v {
                        for k in 0..=dw + r_w {
                            ex[i][j][k][l_new] = ex[i][j][k][l_new].clone()
                                + ew[i][j][k][l_old].clone() * delta.clone();
                        }
                    }
                }
            }
        }
        Self::new(ex)
    }

    pub fn subdivide_u(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();
        let dx = self.coefficients[0][0][0].len();
        let mut left = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        let mut right = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        for j in 0..dv {
            for k in 0..dw {
                for l in 0..dx {
                    let mut beta: Vec<T> = (0..du)
                        .map(|i| self.coefficients[i][j][k][l].clone())
                        .collect();
                    left[0][j][k][l] = beta[0].clone();
                    right[du - 1][j][k][l] = beta[du - 1].clone();
                    for r in 1..du {
                        for i in 0..(du - r) {
                            beta[i] = beta[i].clone() * (EFloat64::one() - t.clone())
                                + beta[i + 1].clone() * t.clone();
                        }
                        left[r][j][k][l] = beta[0].clone();
                        right[du - r - 1][j][k][l] = beta[du - r - 1].clone();
                    }
                }
            }
        }
        (Self::new(left), Self::new(right))
    }

    pub fn subdivide_v(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();
        let dx = self.coefficients[0][0][0].len();
        let mut left = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        let mut right = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        for i in 0..du {
            for k in 0..dw {
                for l in 0..dx {
                    let mut beta: Vec<T> = (0..dv)
                        .map(|j| self.coefficients[i][j][k][l].clone())
                        .collect();
                    left[i][0][k][l] = beta[0].clone();
                    right[i][dv - 1][k][l] = beta[dv - 1].clone();
                    for r in 1..dv {
                        for j in 0..(dv - r) {
                            beta[j] = beta[j].clone() * (EFloat64::one() - t.clone())
                                + beta[j + 1].clone() * t.clone();
                        }
                        left[i][r][k][l] = beta[0].clone();
                        right[i][dv - r - 1][k][l] = beta[dv - r - 1].clone();
                    }
                }
            }
        }
        (Self::new(left), Self::new(right))
    }

    pub fn subdivide_w(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();
        let dx = self.coefficients[0][0][0].len();
        let mut left = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        let mut right = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        for i in 0..du {
            for j in 0..dv {
                for l in 0..dx {
                    let mut beta: Vec<T> = (0..dw)
                        .map(|k| self.coefficients[i][j][k][l].clone())
                        .collect();
                    left[i][j][0][l] = beta[0].clone();
                    right[i][j][dw - 1][l] = beta[dw - 1].clone();
                    for r in 1..dw {
                        for k in 0..(dw - r) {
                            beta[k] = beta[k].clone() * (EFloat64::one() - t.clone())
                                + beta[k + 1].clone() * t.clone();
                        }
                        left[i][j][r][l] = beta[0].clone();
                        right[i][j][dw - r - 1][l] = beta[dw - r - 1].clone();
                    }
                }
            }
        }
        (Self::new(left), Self::new(right))
    }

    pub fn subdivide_x(&self, t: EFloat64) -> (Self, Self) {
        let du = self.coefficients.len();
        let dv = self.coefficients[0].len();
        let dw = self.coefficients[0][0].len();
        let dx = self.coefficients[0][0][0].len();
        let mut left = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        let mut right = vec![vec![vec![vec![T::zero(); dx]; dw]; dv]; du];
        for i in 0..du {
            for j in 0..dv {
                for k in 0..dw {
                    let mut beta: Vec<T> = (0..dx)
                        .map(|l| self.coefficients[i][j][k][l].clone())
                        .collect();
                    left[i][j][k][0] = beta[0].clone();
                    right[i][j][k][dx - 1] = beta[dx - 1].clone();
                    for r in 1..dx {
                        for l in 0..(dx - r) {
                            beta[l] = beta[l].clone() * (EFloat64::one() - t.clone())
                                + beta[l + 1].clone() * t.clone();
                        }
                        left[i][j][k][r] = beta[0].clone();
                        right[i][j][k][dx - r - 1] = beta[dx - r - 1].clone();
                    }
                }
            }
        }
        (Self::new(left), Self::new(right))
    }
}

impl<T> BernsteinHyperVolume<T>
where
    T: Zero + Clone + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn derivative_u(&self) -> BernsteinHyperVolume<T> {
        let du = self.degree_u();
        if du == 0 {
            return BernsteinHyperVolume::new(vec![self.coefficients[0].clone()]);
        }
        let dv = self.degree_v();
        let dw = self.degree_w();
        let dx = self.degree_x();
        let scale = EFloat64::from(du as f64);
        let mut coeffs = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; dv + 1]; du];
        for i in 0..du {
            for j in 0..=dv {
                for k in 0..=dw {
                    for l in 0..=dx {
                        let diff = self.coefficients[i + 1][j][k][l].clone()
                            - self.coefficients[i][j][k][l].clone();
                        coeffs[i][j][k][l] = diff * scale.clone();
                    }
                }
            }
        }
        BernsteinHyperVolume::new(coeffs)
    }

    pub fn derivative_v(&self) -> BernsteinHyperVolume<T> {
        let dv = self.degree_v();
        if dv == 0 {
            let du = self.degree_u();
            let dw = self.degree_w();
            let dx = self.degree_x();
            let mut base = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; 1]; du + 1];
            for i in 0..=du {
                for k in 0..=dw {
                    for l in 0..=dx {
                        base[i][0][k][l] = self.coefficients[i][0][k][l].clone();
                    }
                }
            }
            return BernsteinHyperVolume::new(base);
        }
        let du = self.degree_u();
        let dw = self.degree_w();
        let dx = self.degree_x();
        let scale = EFloat64::from(dv as f64);
        let mut coeffs = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; dv]; du + 1];
        for i in 0..=du {
            for j in 0..dv {
                for k in 0..=dw {
                    for l in 0..=dx {
                        let diff = self.coefficients[i][j + 1][k][l].clone()
                            - self.coefficients[i][j][k][l].clone();
                        coeffs[i][j][k][l] = diff * scale.clone();
                    }
                }
            }
        }
        BernsteinHyperVolume::new(coeffs)
    }

    pub fn derivative_w(&self) -> BernsteinHyperVolume<T> {
        let dw = self.degree_w();
        if dw == 0 {
            let du = self.degree_u();
            let dv = self.degree_v();
            let dx = self.degree_x();
            let mut base = vec![vec![vec![vec![T::zero(); dx + 1]; 1]; dv + 1]; du + 1];
            for i in 0..=du {
                for j in 0..=dv {
                    for l in 0..=dx {
                        base[i][j][0][l] = self.coefficients[i][j][0][l].clone();
                    }
                }
            }
            return BernsteinHyperVolume::new(base);
        }
        let du = self.degree_u();
        let dv = self.degree_v();
        let dx = self.degree_x();
        let scale = EFloat64::from(dw as f64);
        let mut coeffs = vec![vec![vec![vec![T::zero(); dx + 1]; dw]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..dw {
                    for l in 0..=dx {
                        let diff = self.coefficients[i][j][k + 1][l].clone()
                            - self.coefficients[i][j][k][l].clone();
                        coeffs[i][j][k][l] = diff * scale.clone();
                    }
                }
            }
        }
        BernsteinHyperVolume::new(coeffs)
    }

    pub fn derivative_x(&self) -> BernsteinHyperVolume<T> {
        let dx = self.degree_x();
        if dx == 0 {
            let du = self.degree_u();
            let dv = self.degree_v();
            let dw = self.degree_w();
            let mut base = vec![vec![vec![vec![T::zero(); 1]; dw + 1]; dv + 1]; du + 1];
            for i in 0..=du {
                for j in 0..=dv {
                    for k in 0..=dw {
                        base[i][j][k][0] = self.coefficients[i][j][k][0].clone();
                    }
                }
            }
            return BernsteinHyperVolume::new(base);
        }
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        let scale = EFloat64::from(dx as f64);
        let mut coeffs = vec![vec![vec![vec![T::zero(); dx]; dw + 1]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..=dw {
                    for l in 0..dx {
                        let diff = self.coefficients[i][j][k][l + 1].clone()
                            - self.coefficients[i][j][k][l].clone();
                        coeffs[i][j][k][l] = diff * scale.clone();
                    }
                }
            }
        }
        BernsteinHyperVolume::new(coeffs)
    }
}

impl<T> BernsteinHyperVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    fn equalize_degree(
        mut self,
        mut rhs: BernsteinHyperVolume<T>,
    ) -> (BernsteinHyperVolume<T>, BernsteinHyperVolume<T>) {
        let (du1, dv1, dw1, dx1) = (
            self.degree_u(),
            self.degree_v(),
            self.degree_w(),
            self.degree_x(),
        );
        let (du2, dv2, dw2, dx2) = (
            rhs.degree_u(),
            rhs.degree_v(),
            rhs.degree_w(),
            rhs.degree_x(),
        );
        if du1 < du2 {
            self = self.elevate_degree(du2 - du1, 0, 0, 0);
        }
        if dv1 < dv2 {
            self = self.elevate_degree(0, dv2 - dv1, 0, 0);
        }
        if dw1 < dw2 {
            self = self.elevate_degree(0, 0, dw2 - dw1, 0);
        }
        if dx1 < dx2 {
            self = self.elevate_degree(0, 0, 0, dx2 - dx1);
        }
        let (du1, dv1, dw1, dx1) = (
            self.degree_u(),
            self.degree_v(),
            self.degree_w(),
            self.degree_x(),
        );
        if du2 < du1 {
            rhs = rhs.elevate_degree(du1 - du2, 0, 0, 0);
        }
        if dv2 < dv1 {
            rhs = rhs.elevate_degree(0, dv1 - dv2, 0, 0);
        }
        if dw2 < dw1 {
            rhs = rhs.elevate_degree(0, 0, dw1 - dw2, 0);
        }
        if dx2 < dx1 {
            rhs = rhs.elevate_degree(0, 0, 0, dx1 - dx2);
        }
        (self, rhs)
    }
}

impl<T> Add for BernsteinHyperVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let du = lhs.degree_u();
        let dv = lhs.degree_v();
        let dw = lhs.degree_w();
        let dx = lhs.degree_x();
        let mut coeffs = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..=dw {
                    for l in 0..=dx {
                        coeffs[i][j][k][l] = lhs.coefficients[i][j][k][l].clone()
                            + rhs.coefficients[i][j][k][l].clone();
                    }
                }
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Sub for BernsteinHyperVolume<T>
where
    T: Zero + Clone + Sub<Output = T> + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let du = lhs.degree_u();
        let dv = lhs.degree_v();
        let dw = lhs.degree_w();
        let dx = lhs.degree_x();
        let mut coeffs = vec![vec![vec![vec![T::zero(); dx + 1]; dw + 1]; dv + 1]; du + 1];
        for i in 0..=du {
            for j in 0..=dv {
                for k in 0..=dw {
                    for l in 0..=dx {
                        coeffs[i][j][k][l] = lhs.coefficients[i][j][k][l].clone()
                            - rhs.coefficients[i][j][k][l].clone();
                    }
                }
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Mul<EFloat64> for BernsteinHyperVolume<T>
where
    T: Clone + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: EFloat64) -> Self::Output {
        let mut coeffs = self.coefficients.clone();
        for a in &mut coeffs {
            for b in a {
                for c in b {
                    for d in c {
                        *d = d.clone() * rhs;
                    }
                }
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Mul<BernsteinHyperVolume<EFloat64>> for BernsteinHyperVolume<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: BernsteinHyperVolume<EFloat64>) -> Self::Output {
        let (du, dv, dw, dx) = (
            self.degree_u(),
            self.degree_v(),
            self.degree_w(),
            self.degree_x(),
        );
        let (eu, ev, ew, ex) = (
            rhs.degree_u(),
            rhs.degree_v(),
            rhs.degree_w(),
            rhs.degree_x(),
        );
        let mut coeffs =
            vec![vec![vec![vec![T::zero(); dx + ex + 1]; dw + ew + 1]; dv + ev + 1]; du + eu + 1];
        for p in 0..=du + eu {
            let i_min = p.saturating_sub(eu);
            let i_max = du.min(p);
            for q in 0..=dv + ev {
                let j_min = q.saturating_sub(ev);
                let j_max = dv.min(q);
                for r in 0..=dw + ew {
                    let k_min = r.saturating_sub(ew);
                    let k_max = dw.min(r);
                    for s in 0..=dx + ex {
                        let l_min = s.saturating_sub(ex);
                        let l_max = dx.min(s);
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
                                    for l in l_min..=l_max {
                                        let d = s - l;
                                        let fx =
                                            (EFloat64::from(binomial_coefficient(dx, l) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(ex, d) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(dx + ex, s) as f64
                                                ))
                                            .unwrap();
                                        let factor = fu * fv * fw * fx;
                                        acc = acc
                                            + self.coefficients[i][j][k][l].clone()
                                                * (rhs.coefficients[a][b][c][d].clone() * factor);
                                    }
                                }
                            }
                        }
                        coeffs[p][q][r][s] = acc;
                    }
                }
            }
        }
        Self::new(coeffs)
    }
}

impl BernsteinHyperVolume<Point> {
    pub fn dot(&self, rhs: &Self) -> BernsteinHyperVolume<EFloat64> {
        let (du, dv, dw, dx) = (
            self.degree_u(),
            self.degree_v(),
            self.degree_w(),
            self.degree_x(),
        );
        let (eu, ev, ew, ex) = (
            rhs.degree_u(),
            rhs.degree_v(),
            rhs.degree_w(),
            rhs.degree_x(),
        );
        let mut coeffs =
            vec![
                vec![vec![vec![EFloat64::zero(); dx + ex + 1]; dw + ew + 1]; dv + ev + 1];
                du + eu + 1
            ];
        for p in 0..=du + eu {
            let i_min = p.saturating_sub(eu);
            let i_max = du.min(p);
            for q in 0..=dv + ev {
                let j_min = q.saturating_sub(ev);
                let j_max = dv.min(q);
                for r in 0..=dw + ew {
                    let k_min = r.saturating_sub(ew);
                    let k_max = dw.min(r);
                    for s in 0..=dx + ex {
                        let l_min = s.saturating_sub(ex);
                        let l_max = dx.min(s);
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
                                    for l in l_min..=l_max {
                                        let d = s - l;
                                        let fx =
                                            (EFloat64::from(binomial_coefficient(dx, l) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(ex, d) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(dx + ex, s) as f64
                                                ))
                                            .unwrap();
                                        let factor = fu * fv * fw * fx;
                                        acc = acc
                                            + self.coefficients[i][j][k][l]
                                                .dot(rhs.coefficients[a][b][c][d])
                                                * factor;
                                    }
                                }
                            }
                        }
                        coeffs[p][q][r][s] = acc;
                    }
                }
            }
        }
        BernsteinHyperVolume::new(coeffs)
    }

    pub fn cross(&self, rhs: &Self) -> Self {
        let (du, dv, dw, dx) = (
            self.degree_u(),
            self.degree_v(),
            self.degree_w(),
            self.degree_x(),
        );
        let (eu, ev, ew, ex) = (
            rhs.degree_u(),
            rhs.degree_v(),
            rhs.degree_w(),
            rhs.degree_x(),
        );
        let mut coeffs = vec![
            vec![vec![vec![Point::zero(); dx + ex + 1]; dw + ew + 1]; dv + ev + 1];
            du + eu + 1
        ];
        for p in 0..=du + eu {
            let i_min = p.saturating_sub(eu);
            let i_max = du.min(p);
            for q in 0..=dv + ev {
                let j_min = q.saturating_sub(ev);
                let j_max = dv.min(q);
                for r in 0..=dw + ew {
                    let k_min = r.saturating_sub(ew);
                    let k_max = dw.min(r);
                    for s in 0..=dx + ex {
                        let l_min = s.saturating_sub(ex);
                        let l_max = dx.min(s);
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
                                    for l in l_min..=l_max {
                                        let d = s - l;
                                        let fx =
                                            (EFloat64::from(binomial_coefficient(dx, l) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(ex, d) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(dx + ex, s) as f64
                                                ))
                                            .unwrap();
                                        let factor = fu * fv * fw * fx;
                                        acc = acc
                                            + self.coefficients[i][j][k][l]
                                                .cross(rhs.coefficients[a][b][c][d])
                                                * factor;
                                    }
                                }
                            }
                        }
                        coeffs[p][q][r][s] = acc;
                    }
                }
            }
        }
        Self::new(coeffs)
    }
}

impl<T> Div<BernsteinHyperVolume<EFloat64>> for BernsteinHyperVolume<T>
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
    fn div(self, rhs: BernsteinHyperVolume<EFloat64>) -> AlgebraResult<Self> {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        let dx = self.degree_x();
        let eu = rhs.degree_u();
        let ev = rhs.degree_v();
        let ew = rhs.degree_w();
        let ex = rhs.degree_x();
        if du < eu || dv < ev || dw < ew || dx < ex {
            return Err("Division degree mismatch: dividend degree < divisor degree".into());
        }

        let qu = du - eu;
        let qv = dv - ev;
        let qw = dw - ew;
        let qx = dx - ex;
        let mut q = vec![vec![vec![vec![T::zero(); qx + 1]; qw + 1]; qv + 1]; qu + 1];

        let r0000 = rhs.coefficients[0][0][0][0];
        if r0000 == 0.0 {
            return Err("Division by hypervolume with zero at (0,0,0,0) is not supported".into());
        }

        for a in 0..=qu {
            for b in 0..=qv {
                for c in 0..=qw {
                    for d in 0..=qx {
                        let mut s_acc = T::zero();
                        for i in 0..=a {
                            for j in 0..=b {
                                for k in 0..=c {
                                    for l in 0..=d {
                                        if i == a && j == b && k == c && l == d {
                                            continue;
                                        }
                                        if a - i > eu || b - j > ev || c - k > ew || d - l > ex {
                                            continue;
                                        }
                                        let au =
                                            (EFloat64::from(binomial_coefficient(qu, i) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(eu, a - i) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(qu + eu, a) as f64
                                                ))
                                            .unwrap();
                                        let av =
                                            (EFloat64::from(binomial_coefficient(qv, j) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(ev, b - j) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(qv + ev, b) as f64
                                                ))
                                            .unwrap();
                                        let aw =
                                            (EFloat64::from(binomial_coefficient(qw, k) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(ew, c - k) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(qw + ew, c) as f64
                                                ))
                                            .unwrap();
                                        let ax =
                                            (EFloat64::from(binomial_coefficient(qx, l) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(ex, d - l) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(qx + ex, d) as f64
                                                ))
                                            .unwrap();
                                        let alpha = au * av * aw * ax;
                                        s_acc = s_acc
                                            + q[i][j][k][l].clone()
                                                * (rhs.coefficients[a - i][b - j][c - k][d - l]
                                                    * alpha);
                                    }
                                }
                            }
                        }

                        let ak = (EFloat64::from(binomial_coefficient(qu, a) as f64)
                            / EFloat64::from(binomial_coefficient(qu + eu, a) as f64))
                        .unwrap();
                        let al = (EFloat64::from(binomial_coefficient(qv, b) as f64)
                            / EFloat64::from(binomial_coefficient(qv + ev, b) as f64))
                        .unwrap();
                        let am = (EFloat64::from(binomial_coefficient(qw, c) as f64)
                            / EFloat64::from(binomial_coefficient(qw + ew, c) as f64))
                        .unwrap();
                        let an = (EFloat64::from(binomial_coefficient(qx, d) as f64)
                            / EFloat64::from(binomial_coefficient(qx + ex, d) as f64))
                        .unwrap();
                        let denom = r0000 * (ak * al * am * an);

                        let numer = self.coefficients[a][b][c][d].clone() - s_acc;
                        q[a][b][c][d] = (numer / denom)?;
                    }
                }
            }
        }

        let quotient = BernsteinHyperVolume::new(q);
        let recomposed: BernsteinHyperVolume<T> = quotient.clone() * rhs;
        if recomposed != self {
            return Err("Division has a remainder; not exactly divisible".into());
        }
        Ok(quotient)
    }
}

impl BernsteinHyperVolume<EFloat64> {
    pub fn sqrt(&self) -> AlgebraResult<BernsteinHyperVolume<EFloat64>> {
        let du = self.degree_u();
        let dv = self.degree_v();
        let dw = self.degree_w();
        let dx = self.degree_x();
        if du % 2 == 1 || dv % 2 == 1 || dw % 2 == 1 || dx % 2 == 1 {
            return Err(
                "Square root only defined for even degrees (perfect square) hypervolumes".into(),
            );
        }
        let mu = du / 2;
        let mv = dv / 2;
        let mw = dw / 2;
        let mx = dx / 2;

        let c = &self.coefficients;
        if c[0][0][0][0] == 0.0 {
            if du >= 1 && c[1][0][0][0] != 0.0 {
                return Err("Not a perfect square (odd first non-zero in u)".into());
            }
            if dv >= 1 && c[0][1][0][0] != 0.0 {
                return Err("Not a perfect square (odd first non-zero in v)".into());
            }
            if dw >= 1 && c[0][0][1][0] != 0.0 {
                return Err("Not a perfect square (odd first non-zero in w)".into());
            }
            if dx >= 1 && c[0][0][0][1] != 0.0 {
                return Err("Not a perfect square (odd first non-zero in x)".into());
            }
        }

        let d0000 = c[0][0][0][0].sqrt()?;
        let mut d = vec![vec![vec![vec![EFloat64::zero(); mx + 1]; mw + 1]; mv + 1]; mu + 1];
        d[0][0][0][0] = d0000;

        for s in 0..=(mu + mv + mw + mx) {
            for a in 0..=s.min(mu) {
                for b in 0..=(s - a).min(mv) {
                    for cidx in 0..=((s - a - b).min(mw)) {
                        let didx = s - a - b - cidx;
                        if didx > mx {
                            continue;
                        }
                        if a == 0 && b == 0 && cidx == 0 && didx == 0 {
                            continue;
                        }
                        let mut s_acc = EFloat64::zero();
                        for i in 0..=a {
                            for j in 0..=b {
                                for k in 0..=cidx {
                                    for l in 0..=didx {
                                        if i == a && j == b && k == cidx && l == didx {
                                            continue;
                                        }
                                        let fu =
                                            (EFloat64::from(binomial_coefficient(mu, i) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(mu, a - i) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(2 * mu, a) as f64
                                                ))
                                            .unwrap();
                                        let fv =
                                            (EFloat64::from(binomial_coefficient(mv, j) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(mv, b - j) as f64
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(2 * mv, b) as f64
                                                ))
                                            .unwrap();
                                        let fw =
                                            (EFloat64::from(binomial_coefficient(mw, k) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(mw, cidx - k) as f64,
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(2 * mw, cidx) as f64,
                                                ))
                                            .unwrap();
                                        let fx =
                                            (EFloat64::from(binomial_coefficient(mx, l) as f64)
                                                * EFloat64::from(
                                                    binomial_coefficient(mx, didx - l) as f64,
                                                )
                                                / EFloat64::from(
                                                    binomial_coefficient(2 * mx, didx) as f64,
                                                ))
                                            .unwrap();
                                        s_acc = s_acc
                                            + (d[i][j][k][l]
                                                * d[a - i][b - j][cidx - k][didx - l]
                                                * (fu * fv * fw * fx));
                                    }
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
                        let ax = (EFloat64::from(binomial_coefficient(mx, didx) as f64)
                            / EFloat64::from(binomial_coefficient(2 * mx, didx) as f64))
                        .unwrap();
                        let denom = EFloat64::from(2.0) * au * av * aw * ax * d[0][0][0][0];
                        let cabcd = c[a][b][cidx][didx];
                        d[a][b][cidx][didx] = ((cabcd - s_acc) / denom)?;
                    }
                }
            }
        }

        let q = BernsteinHyperVolume::new(d);
        let back = q.clone() * q.clone();
        if back != *self {
            return Err("Verification of hypervolume square root failed".into());
        }
        Ok(q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_hyper() -> BernsteinHyperVolume<Point> {
        // 2x2x2x2 grid (degree 1,1,1,1) simple affine layout
        let mut coeffs = vec![vec![vec![vec![Point::zero(); 2]; 2]; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        coeffs[i][j][k][l] = Point::from_f64(
                            i as f64,
                            j as f64 + k as f64 * 0.1,
                            k as f64 + l as f64 * 0.2,
                        );
                    }
                }
            }
        }
        BernsteinHyperVolume::new(coeffs)
    }

    #[test]
    fn test_eval_center() {
        let h = sample_hyper();
        let p = h.eval(
            EFloat64::from(0.5),
            EFloat64::from(0.5),
            EFloat64::from(0.5),
            EFloat64::from(0.5),
        );
        assert!(p != Point::zero());
    }

    #[test]
    fn test_elevate_and_eval_invariance() {
        let h = sample_hyper();
        let h2 = h.elevate_degree(1, 1, 1, 1);
        for iu in 0..=2 {
            for iv in 0..=2 {
                for iw in 0..=2 {
                    for ix in 0..=2 {
                        let u = EFloat64::from(iu as f64 / 2.0);
                        let v = EFloat64::from(iv as f64 / 2.0);
                        let w = EFloat64::from(iw as f64 / 2.0);
                        let x = EFloat64::from(ix as f64 / 2.0);
                        assert_eq!(
                            h.eval(u.clone(), v.clone(), w.clone(), x.clone()),
                            h2.eval(u, v, w, x)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_subdivide_all() {
        let h = sample_hyper();
        let t = EFloat64::from(0.5);
        let (lu, ru) = h.subdivide_u(t.clone());
        let (lv, rv) = h.subdivide_v(t.clone());
        let (lw, rw) = h.subdivide_w(t.clone());
        let (lx, rx) = h.subdivide_x(t.clone());

        for iu in 0..=2 {
            for iv in 0..=2 {
                for iw in 0..=2 {
                    for ix in 0..=2 {
                        let u = EFloat64::from(iu as f64 / 2.0);
                        let v = EFloat64::from(iv as f64 / 2.0);
                        let w = EFloat64::from(iw as f64 / 2.0);
                        let x = EFloat64::from(ix as f64 / 2.0);
                        // u-left mapping: (u*t, v, w, x)
                        assert_eq!(
                            h.eval(u.clone() * t.clone(), v.clone(), w.clone(), x.clone()),
                            lu.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // u-right mapping: (t + u*(1-t), v, w, x)
                        let one_minus_t = EFloat64::one() - t.clone();
                        let mu = t.clone() + u.clone() * one_minus_t.clone();
                        assert_eq!(
                            h.eval(mu, v.clone(), w.clone(), x.clone()),
                            ru.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // v-left
                        assert_eq!(
                            h.eval(u.clone(), v.clone() * t.clone(), w.clone(), x.clone()),
                            lv.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // v-right
                        let mv = t.clone() + v.clone() * one_minus_t.clone();
                        assert_eq!(
                            h.eval(u.clone(), mv, w.clone(), x.clone()),
                            rv.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // w-left
                        assert_eq!(
                            h.eval(u.clone(), v.clone(), w.clone() * t.clone(), x.clone()),
                            lw.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // w-right
                        let mw = t.clone() + w.clone() * one_minus_t.clone();
                        assert_eq!(
                            h.eval(u.clone(), v.clone(), mw, x.clone()),
                            rw.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // x-left
                        assert_eq!(
                            h.eval(u.clone(), v.clone(), w.clone(), x.clone() * t.clone()),
                            lx.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                        // x-right
                        let mx = t.clone() + x.clone() * one_minus_t.clone();
                        assert_eq!(
                            h.eval(u.clone(), v.clone(), w.clone(), mx),
                            rx.eval(u.clone(), v.clone(), w.clone(), x.clone())
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_derivatives() {
        let h = sample_hyper();
        let hu = h.derivative_u();
        let hv = h.derivative_v();
        let hw = h.derivative_w();
        let hx = h.derivative_x();
        let u = EFloat64::from(0.3);
        let v = EFloat64::from(0.4);
        let w = EFloat64::from(0.5);
        let x = EFloat64::from(0.6);
        let _ = hu.eval(u.clone(), v.clone(), w.clone(), x.clone());
        let _ = hv.eval(u.clone(), v.clone(), w.clone(), x.clone());
        let _ = hw.eval(u.clone(), v.clone(), w.clone(), x.clone());
        let _ = hx.eval(u, v, w, x);
    }

    #[test]
    fn test_add_sub_scalar_mul() {
        let h1 = sample_hyper();
        let h2 = sample_hyper();
        let sum = h1.clone() + h2.clone();
        let diff = h1.clone() - h2.clone();
        let sm = h1.clone() * EFloat64::from(2.5);
        let u = EFloat64::from(0.25);
        let v = EFloat64::from(0.5);
        let w = EFloat64::from(0.75);
        let x = EFloat64::from(0.33);
        assert_eq!(
            sum.eval(u.clone(), v.clone(), w.clone(), x.clone()),
            h1.eval(u.clone(), v.clone(), w.clone(), x.clone())
                + h2.eval(u.clone(), v.clone(), w.clone(), x.clone())
        );
        assert_eq!(
            diff.eval(u.clone(), v.clone(), w.clone(), x.clone()),
            h1.eval(u.clone(), v.clone(), w.clone(), x.clone())
                - h2.eval(u.clone(), v.clone(), w.clone(), x.clone())
        );
        assert_eq!(
            sm.eval(u.clone(), v.clone(), w.clone(), x.clone()),
            h1.eval(u.clone(), v.clone(), w.clone(), x.clone()) * EFloat64::from(2.5)
        );
    }

    #[test]
    fn test_dot_cross_mul_div_sqrt() -> AlgebraResult<()> {
        let h = sample_hyper();
        let h2 = sample_hyper();
        let dot = h.dot(&h2);
        let cross = h.cross(&h2);

        // numeric mul/div roundtrip
        let q = BernsteinHyperVolume::new(vec![vec![vec![vec![
            EFloat64::from(0.5),
            EFloat64::from(0.5),
        ]]]]); // degree (0,0,0,1)
        let p = h.clone() * q.clone();
        let r = (p.clone() / q.clone())?;
        assert_eq!(r, h);

        // dot/cross equality at sample point
        let u = EFloat64::from(0.2);
        let v = EFloat64::from(0.4);
        let w = EFloat64::from(0.6);
        let x = EFloat64::from(0.8);
        let lhs_dot = h
            .eval(u.clone(), v.clone(), w.clone(), x.clone())
            .dot(h2.eval(u.clone(), v.clone(), w.clone(), x.clone()));
        let rhs_dot = dot.eval(u.clone(), v.clone(), w.clone(), x.clone());
        assert_eq!(lhs_dot, rhs_dot);
        let lhs_cross = h
            .eval(u.clone(), v.clone(), w.clone(), x.clone())
            .cross(h2.eval(u.clone(), v.clone(), w.clone(), x.clone()));
        let rhs_cross = cross.eval(u.clone(), v.clone(), w.clone(), x.clone());
        assert_eq!(lhs_cross, rhs_cross);

        // sqrt: construct q (1,1,1,1), p = q*q then sqrt
        let q2 = BernsteinHyperVolume::new(vec![
            vec![
                vec![
                    vec![EFloat64::from(0.5), EFloat64::from(1.0)],
                    vec![EFloat64::from(1.0), EFloat64::from(2.0)],
                ],
                vec![
                    vec![EFloat64::from(1.0), EFloat64::from(2.0)],
                    vec![EFloat64::from(2.0), EFloat64::from(3.0)],
                ],
            ],
            vec![
                vec![
                    vec![EFloat64::from(1.0), EFloat64::from(1.5)],
                    vec![EFloat64::from(2.0), EFloat64::from(2.5)],
                ],
                vec![
                    vec![EFloat64::from(2.0), EFloat64::from(2.5)],
                    vec![EFloat64::from(3.0), EFloat64::from(3.5)],
                ],
            ],
        ]);
        let p2 = q2.clone() * q2.clone();
        let srt = p2.sqrt()?;
        let back = srt.clone() * srt.clone();
        assert_eq!(back, p2);

        Ok(())
    }
}
