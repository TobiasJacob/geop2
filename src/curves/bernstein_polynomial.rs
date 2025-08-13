use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

use crate::binomial_coefficient;
use crate::zero::Zero;
use crate::{
    algebra_error::AlgebraResult,
    curves::curve_like::CurveLike,
    primitives::{convex_hull::ConvexHull, efloat::EFloat64, point::Point},
};

// Represents a polynomial in the form of a_{0} B_{0,n}
#[derive(Debug, Clone)]
pub struct BernsteinPolynomial<T> {
    pub coefficients: Vec<T>,
}

impl<T> BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn bernstein_basis(i: usize, n: usize, up_vector: T, direction_vector: T) -> Self {
        let mut coefficients = vec![T::zero(); n + 1];
        for j in 0..=n {
            coefficients[j] = direction_vector.clone() * EFloat64::from(j as f64 / n as f64);
        }
        coefficients[i] = coefficients[i].clone() + up_vector;
        Self::new(coefficients)
    }
}

impl<T> BernsteinPolynomial<T> {
    pub fn new(coefficients: Vec<T>) -> Self {
        Self { coefficients }
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }
}

impl<T> Display for BernsteinPolynomial<T>
where
    T: Display + Zero + PartialEq,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        let n = self.degree();
        for (i, coeff) in self.coefficients.iter().enumerate() {
            if *coeff != T::zero() {
                if !first {
                    write!(f, " + ")?;
                }
                write!(f, "{} B_{{{},{}}}(t)", coeff, i, n)?;
                first = false;
            }
        }
        Ok(())
    }
}

impl<T> BernsteinPolynomial<T>
where
    T: Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn eval(&self, t: EFloat64) -> T {
        let mut beta = self.coefficients.clone();
        let n = beta.len();
        for j in 1..n {
            for k in 0..n - j {
                beta[k] = beta[k].clone() * (EFloat64::one() - t.clone())
                    + beta[k + 1].clone() * t.clone();
            }
        }
        beta[0].clone()
    }
}

impl<T> BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    //$$ c_i^{n+r} = \sum_{j = max(0, i - r)}^{min(n, i)} \frac{\binom{r}{i - j} \binom{n}{j}}{\binom{n + r}{i}} c_i^n $$
    pub fn elevate_degree(&self, r: usize) -> Self {
        let n = self.degree();
        let mut new_coeffs = vec![T::zero(); n + r + 1];

        for i in 0..=n + r {
            for j in i.saturating_sub(r)..=n.min(i) {
                let factor = (EFloat64::from(
                    (binomial_coefficient(r, i - j) * binomial_coefficient(n, j)) as f64,
                ) / EFloat64::from(binomial_coefficient(n + r, i) as f64))
                .unwrap();
                new_coeffs[i] = new_coeffs[i].clone() + (self.coefficients[j].clone() * factor);
            }
        }

        Self::new(new_coeffs)
    }

    pub fn equalize_degree(
        self,
        rhs: BernsteinPolynomial<T>,
    ) -> (BernsteinPolynomial<T>, BernsteinPolynomial<T>) {
        let n = self.degree();
        let m = rhs.degree();
        if n == m {
            return (self, rhs);
        }
        if n < m {
            let lhs = self.elevate_degree(m - n);
            (lhs, rhs)
        } else {
            let rhs = rhs.elevate_degree(n - m);
            (self, rhs)
        }
    }

    // Use de Casteljau's algorithm to subdivide the polynomial
    pub fn subdivide(&self, t: EFloat64) -> (BernsteinPolynomial<T>, BernsteinPolynomial<T>) {
        let mut beta = self.coefficients.clone();
        let n = beta.len();
        let mut left = vec![T::zero(); n];
        let mut right = vec![T::zero(); n];

        left[0] = beta[0].clone();
        right[n - 1] = beta[n - 1].clone();
        for j in 1..n {
            for k in 0..n - j {
                beta[k] = beta[k].clone() * (EFloat64::one() - t.clone())
                    + beta[k + 1].clone() * t.clone();
            }
            left[j] = beta[0].clone();
            right[n - j - 1] = beta[n - j - 1].clone();
        }

        (Self::new(left), Self::new(right))
    }
}

impl<T> BernsteinPolynomial<T>
where
    T: Zero
        + Clone
        + Sub<Output = T>
        + Mul<EFloat64, Output = T>
        + Div<EFloat64, Output = AlgebraResult<T>>
        + PartialEq,
{
    // Try to reduce the degree by exactly one. If not possible, return an error.
    fn reduce_degree_once(&self) -> AlgebraResult<Self> {
        let n = self.degree();
        if n == 0 {
            return Err("Cannot reduce degree of a constant polynomial".into());
        }

        // q will have degree n-1
        let mut q: Vec<T> = Vec::with_capacity(n);

        // Boundary condition from elevation inverse: q_0 = a_0
        q.push(self.coefficients[0].clone());

        let n_ef = EFloat64::from(n as f64);
        for i in 1..n {
            // q_i = (n * a_i - i * q_{i-1}) / (n - i)
            let numer = self.coefficients[i].clone() * n_ef.clone()
                - q[i - 1].clone() * EFloat64::from(i as f64);
            let denom = EFloat64::from((n - i) as f64);
            q.push((numer / denom)?);
        }

        // Consistency check with the second boundary condition: q_{n-1} must equal a_n
        if q[n - 1] != self.coefficients[n] {
            return Err("Degree reduction not possible exactly for given coefficients".into());
        }

        Ok(Self::new(q))
    }

    // Reduce degree by r steps. If any step fails, return error.
    pub fn reduce_degree(&self, r: usize) -> AlgebraResult<Self> {
        if r == 0 {
            return Ok(self.clone());
        }
        if r > self.degree() {
            return Err("Reduction degree exceeds polynomial degree".into());
        }

        let mut current = self.clone();
        for _ in 0..r {
            current = current.reduce_degree_once()?;
        }
        Ok(current)
    }
}

impl<T> BernsteinPolynomial<T>
where
    T: Zero + Clone + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    pub fn derivative(&self) -> BernsteinPolynomial<T> {
        let n = self.degree();
        if n == 0 {
            return BernsteinPolynomial::new(vec![T::zero()]);
        }

        let scale = EFloat64::from(n as f64);
        let mut deriv_coeffs = Vec::with_capacity(n);
        for i in 0..n {
            let diff = self.coefficients[i + 1].clone() - self.coefficients[i].clone();
            deriv_coeffs.push(diff * scale);
        }
        BernsteinPolynomial::new(deriv_coeffs)
    }
}
impl<T> PartialEq for BernsteinPolynomial<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.coefficients == other.coefficients
    }
}

impl<T> Add for BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let coefficients = lhs
            .coefficients
            .into_iter()
            .zip(rhs.coefficients.into_iter())
            .map(|(a, b)| a + b)
            .collect();
        return Self::new(coefficients);
    }
}

impl<T> Sub for BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = self.equalize_degree(rhs);
        let coefficients = lhs
            .coefficients
            .into_iter()
            .zip(rhs.coefficients.into_iter())
            .map(|(a, b)| a - b)
            .collect();
        return Self::new(coefficients);
    }
}

// Scalar multiplication of a Bernstein polynomial by EFloat64
impl<T> Mul<EFloat64> for BernsteinPolynomial<T>
where
    T: Clone + Mul<EFloat64, Output = T>,
{
    type Output = Self;

    fn mul(self, rhs: EFloat64) -> Self::Output {
        let coefficients = self.coefficients.into_iter().map(|c| c * rhs).collect();
        Self::new(coefficients)
    }
}

// Multiplication of two Bernstein polynomials with numeric coefficients
impl<T> Mul<BernsteinPolynomial<EFloat64>> for BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;

    fn mul(self, rhs: BernsteinPolynomial<EFloat64>) -> Self::Output {
        let n = self.degree();
        let m = rhs.degree();
        let mut coefficients = vec![T::zero(); n + m + 1];

        for k in 0..=n + m {
            let i_min = k.saturating_sub(m);
            let i_max = n.min(k);
            let mut acc = T::zero();
            for i in i_min..=i_max {
                let j = k - i;
                let factor = (EFloat64::from(binomial_coefficient(n, i) as f64)
                    * EFloat64::from(binomial_coefficient(m, j) as f64)
                    / EFloat64::from(binomial_coefficient(n + m, k) as f64))
                .unwrap();
                acc = acc + self.coefficients[i].clone() * rhs.coefficients[j].clone() * factor;
            }
            coefficients[k] = acc;
        }

        Self::new(coefficients)
    }
}

// Exact division by a numeric Bernstein polynomial. Errors if there is a remainder.
impl<T> Div<BernsteinPolynomial<EFloat64>> for BernsteinPolynomial<T>
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

    fn div(self, rhs: BernsteinPolynomial<EFloat64>) -> AlgebraResult<Self> {
        let p_deg = self.degree();
        let r_deg = rhs.degree();
        if p_deg < r_deg {
            return Err("Division degree mismatch: dividend degree < divisor degree".into());
        }

        let q_deg = p_deg - r_deg;
        let mut q = vec![T::zero(); q_deg + 1];

        // Require r_0 to be non-zero for stable forward substitution
        let r0 = rhs.coefficients[0];
        if r0 == 0.0 {
            return Err("Division by polynomial with zero value at t=0 is not supported".into());
        }

        // Forward substitution to solve for q_k
        for k in 0..=q_deg {
            // Compute S = sum_{i=0..k-1} q_i * r_{k-i} * alpha(i,k-i)
            let mut s = T::zero();
            for i in 0..k {
                let j = k - i;
                if j > r_deg {
                    continue;
                }
                let factor = (EFloat64::from(binomial_coefficient(q_deg, i) as f64)
                    * EFloat64::from(binomial_coefficient(r_deg, j) as f64)
                    / EFloat64::from(binomial_coefficient(q_deg + r_deg, k) as f64))
                .unwrap();
                s = s + q[i].clone() * (rhs.coefficients[j] * factor);
            }

            // denom = r_0 * alpha(k,0)
            let alpha_k0 = (EFloat64::from(binomial_coefficient(q_deg, k) as f64)
                / EFloat64::from(binomial_coefficient(q_deg + r_deg, k) as f64))
            .unwrap();
            let denom = r0 * alpha_k0;

            // q_k = (p_k - S) / denom
            let numer = self.coefficients[k].clone() - s;
            q[k] = (numer / denom)?;
        }

        let quotient = BernsteinPolynomial::new(q);

        // Verify exactness: quotient * rhs must equal self exactly
        let recomposed: BernsteinPolynomial<T> = quotient.clone() * rhs;
        if recomposed != self {
            return Err("Division has a remainder; not exactly divisible".into());
        }

        Ok(quotient)
    }
}

// Dot product of two Bernstein polynomials with numeric coefficients
impl BernsteinPolynomial<Point> {
    pub fn dot(&self, rhs: &Self) -> BernsteinPolynomial<EFloat64> {
        let n = self.degree();
        let m = rhs.degree();
        let mut coefficients = vec![EFloat64::zero(); n + m + 1];

        for k in 0..=n + m {
            let i_min = k.saturating_sub(m);
            let i_max = n.min(k);
            let mut acc = EFloat64::zero();
            for i in i_min..=i_max {
                let j = k - i;
                let factor = (EFloat64::from(binomial_coefficient(n, i) as f64)
                    * EFloat64::from(binomial_coefficient(m, j) as f64)
                    / EFloat64::from(binomial_coefficient(n + m, k) as f64))
                .unwrap();
                acc = acc + self.coefficients[i].dot(rhs.coefficients[j]) * factor;
            }
            coefficients[k] = acc;
        }

        BernsteinPolynomial::<EFloat64>::new(coefficients)
    }

    pub fn cross(&self, rhs: &Self) -> Self {
        let n = self.degree();
        let m = rhs.degree();
        let mut coefficients = vec![Point::zero(); n + m + 1];

        for k in 0..=n + m {
            let i_min = k.saturating_sub(m);
            let i_max = n.min(k);
            let mut acc = Point::zero();
            for i in i_min..=i_max {
                let j = k - i;
                let factor = (EFloat64::from(binomial_coefficient(n, i) as f64)
                    * EFloat64::from(binomial_coefficient(m, j) as f64)
                    / EFloat64::from(binomial_coefficient(n + m, k) as f64))
                .unwrap();
                acc = acc + self.coefficients[i].cross(rhs.coefficients[j]) * factor;
            }
            coefficients[k] = acc;
        }

        Self::new(coefficients)
    }
}

impl BernsteinPolynomial<EFloat64> {
    pub fn sqrt(&self) -> AlgebraResult<BernsteinPolynomial<EFloat64>> {
        let n = self.degree();
        if n % 2 == 1 {
            return Err(
                "Square root only defined for even degree (perfect square) polynomials".into(),
            );
        }
        let m = n / 2;

        // Coefficients of p
        let c = &self.coefficients;

        // Handle leading zero(s): if c[0] == 0, then c[1] must be 0, etc.
        // Simple early contradiction check:
        if c[0] == 0.0 {
            if n >= 1 && c[1] != 0.0 {
                return Err("Polynomial is not a perfect square in Bernstein basis (first non-zero coefficient index is odd)".into());
            }
            // More robust handling of multiple leading zeros could be added here.
        }

        // d_0 = sqrt(c_0) (choose the non-negative branch).
        let d0 = c[0].sqrt()?;
        // If you want the negative branch, flip the sign of all d later.

        // Prepare d[0..=m]
        let mut d = vec![EFloat64::zero(); m + 1];
        d[0] = d0;

        // Solve for d_k for k = 1..m
        for k in 1..=m {
            // Accumulate S = sum_{i=1}^{k-1} [ (C(m,i) C(m,k-i) / C(2m,k)) * d_i d_{k-i} ]
            let mut s = EFloat64::zero();
            for i in 1..k {
                let j = k - i;
                let num = binomial_coefficient(m, i) * binomial_coefficient(m, j);
                let den = binomial_coefficient(2 * m, k);
                let factor = (EFloat64::from(num as f64) / EFloat64::from(den as f64))?;
                s = s + (d[i] * d[j] * factor);
            }

            // Coefficient for the linear term in d_k:
            // 2 * (C(m,k)/C(2m,k)) * d_0
            let a = (EFloat64::from(binomial_coefficient(m, k) as f64)
                / EFloat64::from(binomial_coefficient(2 * m, k) as f64))?;
            let denom = EFloat64::from(2.0) * a * d[0];

            // d_k = (c_k - S) / denom
            let numer = c[k] - s;
            d[k] = (numer / denom)?;
        }

        let q = BernsteinPolynomial::new(d);

        // Consistency check: q*q must reproduce self exactly
        let check = q.clone() * q.clone();
        if check != *self {
            return Err("Verification of square root failed".into());
        }

        Ok(q)
    }
}

impl CurveLike for BernsteinPolynomial<Point> {
    fn span(&self) -> (EFloat64, EFloat64) {
        (EFloat64::zero(), EFloat64::one())
    }

    fn eval(&self, t: EFloat64) -> Point {
        self.eval(t)
    }

    fn subdivide(&self, t: EFloat64) -> AlgebraResult<(Self, Self)> {
        Ok(self.subdivide(t))
    }

    fn split(&self) -> AlgebraResult<(Self, Self)> {
        Ok(self.subdivide(EFloat64::from(0.5)))
    }

    fn control_polygon_hull(&self) -> AlgebraResult<ConvexHull> {
        ConvexHull::try_new(self.coefficients.clone())
    }
}

#[cfg(test)]
mod tests {
    use crate::primitives::{color::Color10, line::Line, primitive_scene::PrimitiveScene};

    use super::*;
    fn test_bernstein_polynomial() -> BernsteinPolynomial<Point> {
        let coeffs = vec![
            // Point::unit_z() * EFloat64::from(1.0),
            // Point::unit_z() * EFloat64::from(2.0),
            // Point::unit_z() * EFloat64::from(1.0),
            // Point::unit_z() * EFloat64::from(5.0),
            // Point::unit_z() * EFloat64::from(3.0),
            Point::from_f64(0.0, 1.0, 0.0),
            Point::from_f64(1.0, 2.0, 0.0),
            Point::from_f64(2.0, 1.0, 4.0),
            Point::from_f64(3.0, 5.0, 0.0),
            Point::from_f64(4.0, 3.0, 0.0),
        ];
        let bernstein = BernsteinPolynomial::new(coeffs.clone());
        bernstein
    }

    #[test]
    fn test_bernstein_with_efloat() {
        let coeffs = vec![EFloat64::from(1.0), EFloat64::from(2.0)];
        let bernstein = BernsteinPolynomial::new(coeffs.clone());
        bernstein.eval(EFloat64::from(0.5));
        bernstein.derivative();
        bernstein.elevate_degree(1);
        bernstein.subdivide(EFloat64::from(0.5));
    }

    #[test]
    fn test_bernstein_eval() {
        let bernstein = test_bernstein_polynomial();

        let t = EFloat64::from(0.5);
        let eval = bernstein.eval(t);
        println!("Bernstein Polynomial: {}", &bernstein);
        println!("Bernstein Polynomial at {}: {}", t, eval);
    }

    #[test]
    fn test_bernstein_elevate_degree() {
        let bernstein = test_bernstein_polynomial();

        let r = 2;
        let elevated_bernstein = bernstein.elevate_degree(r);

        println!("Bernstein Polynomial: {}", &bernstein);
        println!("Elevated Bernstein Polynomial: {}", &elevated_bernstein);

        for t in 0..=10 {
            let t = EFloat64::from(t as f64 / 10.0);
            assert_eq!(bernstein.eval(t), elevated_bernstein.eval(t), "t = {}", t);
        }
    }

    #[test]
    fn test_bernstein_elevate_degree2() {
        let coeffs = vec![
            Point::unit_z() * EFloat64::from(1.0),
            Point::unit_z() * EFloat64::from(2.0),
        ];
        let bernstein = BernsteinPolynomial::new(coeffs.clone());

        println!("Bernstein Polynomial: {}", &bernstein);
        println!(
            "Elevated Bernstein Polynomial: {}",
            &bernstein.elevate_degree(1)
        );
        println!(
            "Elevated Bernstein Polynomial 2: {}",
            &bernstein.elevate_degree(2)
        );
    }

    #[test]
    fn test_bernstein_subdivide() {
        let bernstein = test_bernstein_polynomial();

        let t = EFloat64::from(0.5);
        let (left, right) = bernstein.subdivide(t);

        println!("Bernstein Polynomial: {}", &bernstein);
        println!("Left Bernstein Polynomial: {}", &left);
        println!("Right Bernstein Polynomial: {}", &right);

        for t in 0..=10 {
            let t = EFloat64::from(t as f64 / 10.0);
            assert_eq!(
                bernstein.eval((t / EFloat64::two()).unwrap()),
                left.eval(t),
                "t = {}",
                t
            );
        }

        for t in 0..=10 {
            let t = EFloat64::from(t as f64 / 10.0);
            assert_eq!(
                bernstein.eval(((EFloat64::one() + t) / EFloat64::two()).unwrap()),
                right.eval(t),
                "t = {}",
                t
            );
        }
    }

    fn test_bernstein_polynomial2() -> BernsteinPolynomial<Point> {
        let coeffs = vec![
            Point::from_f64(-1.0, 0.0, 2.0),
            Point::from_f64(0.5, -2.0, 1.0),
            Point::from_f64(1.5, 1.0, -1.0),
            Point::from_f64(0.0, 3.0, 0.0),
        ];
        BernsteinPolynomial::new(coeffs)
    }

    #[test]
    fn test_point_addition() {
        let p = test_bernstein_polynomial();
        let q = test_bernstein_polynomial2();
        let sum = p.clone() + q.clone();

        for i in 0..=10 {
            let t = EFloat64::from(i as f64 / 10.0);
            let lhs = sum.eval(t);
            let rhs = p.eval(t) + q.eval(t);
            assert_eq!(lhs, rhs, "t = {}", t);
        }
    }

    #[test]
    fn test_point_subtraction() {
        let p = test_bernstein_polynomial();
        let q = test_bernstein_polynomial2();
        let diff = p.clone() - q.clone();

        for i in 0..=10 {
            let t = EFloat64::from(i as f64 / 10.0);
            let lhs = diff.eval(t);
            let rhs = p.eval(t) - q.eval(t);
            assert_eq!(lhs, rhs, "t = {}", t);
        }
    }

    #[test]
    fn test_point_scalar_multiplication() {
        let p = test_bernstein_polynomial();
        let s = EFloat64::from(2.5);
        let prod = p.clone() * s;

        for i in 0..=10 {
            let t = EFloat64::from(i as f64 / 10.0);
            let lhs = prod.eval(t);
            let rhs = p.eval(t) * s;
            assert_eq!(lhs, rhs, "t = {}", t);
        }
    }

    #[test]
    fn test_point_dot_product() {
        let p = test_bernstein_polynomial();
        let q = test_bernstein_polynomial2();

        for i in 0..=10 {
            let t = EFloat64::from(i as f64 / 10.0);
            let a = p.eval(t);
            let b = q.eval(t);
            let lhs = a.dot(b);
            let rhs = a.dot(b);
            assert_eq!(lhs, rhs, "t = {}", t);
        }
    }

    #[test]
    fn test_bernstein_sqrt_perfect_square() {
        let q = BernsteinPolynomial::new(vec![
            EFloat64::from(0.5),
            EFloat64::from(1.0),
            EFloat64::from(2.0),
        ]);
        let p = q.clone() * q.clone();
        let s = p.sqrt().expect("sqrt should succeed for perfect square");
        assert_eq!(s.degree(), q.degree());
        let recomposed = s.clone() * s.clone();
        assert!(
            recomposed
                .coefficients
                .iter()
                .zip(p.coefficients.iter())
                .all(|(a, b)| *a == *b)
        );
    }

    #[test]
    fn test_bernstein_sqrt_error_non_square() {
        // Odd degree polynomial cannot be a perfect square in this sense
        let p = BernsteinPolynomial::new(vec![
            EFloat64::from(1.0),
            EFloat64::from(2.0),
            EFloat64::from(3.0),
            EFloat64::from(4.0),
        ]);
        assert!(p.sqrt().is_err());
    }

    #[test]
    fn test_bernstein_derivative() -> AlgebraResult<()> {
        let bernstein = test_bernstein_polynomial();
        let tangent = bernstein.derivative();
        let second_derivative = tangent.derivative();
        println!("Bernstein Polynomial: {}", &bernstein);
        println!("Tangent Bernstein Polynomial: {}", &tangent);
        println!(
            "Second Derivative Bernstein Polynomial: {}",
            &second_derivative
        );

        let mut scene = PrimitiveScene::new();
        scene.add_curve(&bernstein, Color10::Blue)?;

        for i in 0..=10 {
            let t = EFloat64::from(i as f64 / 10.0);
            let eval = bernstein.eval(t);
            let tangent_eval = tangent.eval(t) * EFloat64::from(0.1);
            let line = Line::try_new(eval, tangent_eval + eval)?;
            scene.add_line(line, Color10::Red);
            let line = Line::try_new(
                eval,
                second_derivative.eval(t) * EFloat64::from(0.01) + eval,
            )?;
            scene.add_line(line, Color10::Green);
        }

        scene.save_to_file("test_outputs/bernstein_derivative.html")?;
        Ok(())
    }

    #[test]
    fn test_bernstein_polynomial_cross_product() {
        let p = test_bernstein_polynomial();
        let q = test_bernstein_polynomial2();
        let cross = p.cross(&q);
        println!("Cross product: {}", &cross);

        // sample at ten points and check if the cross product is perpendicular to the original polynomials
        for i in 0..=10 {
            let t = EFloat64::from(i as f64 / 10.0);
            let a = p.eval(t);
            let b = q.eval(t);
            let cross_eval = cross.eval(t);
            let dot_product = a.dot(cross_eval);
            assert!(dot_product == 0.0, "t = {}", t);
            let dot_product = b.dot(cross_eval);
            assert!(dot_product == 0.0, "t = {}", t);
        }
    }

    #[test]
    fn test_bernstein_polynomial_division() -> AlgebraResult<()> {
        // Construct two polynomials p = q * r, then check p / q == r and p / r == q

        // Let q(t) = 1 + 2t + t^2 (degree 2)
        let q = BernsteinPolynomial::new(vec![
            EFloat64::from(1.0),
            EFloat64::from(2.0),
            EFloat64::from(1.0),
        ]);
        // Let r(t) = 2 + 0t + 3t^2 (degree 2)
        let r = BernsteinPolynomial::new(vec![
            EFloat64::from(2.0),
            EFloat64::from(0.0),
            EFloat64::from(3.0),
        ]);
        // Compute p = q * r (degree 4)
        let p = q.clone() * r.clone();

        // Division: p / q == r
        let r2 = (p.clone() / q.clone())?;
        assert_eq!(r2, r);

        // Division: p / r == q
        let q2 = (p.clone() / r.clone())?;
        assert_eq!(q2, q);

        // Division by a polynomial of higher degree should error
        let result = q.clone() / p.clone();
        assert!(result.is_err());

        // Division by a polynomial with zero at t=0 should error
        let zero_at_0 = BernsteinPolynomial::new(vec![
            EFloat64::from(0.0),
            EFloat64::from(1.0),
            EFloat64::from(2.0),
        ]);
        let result = p.clone() / zero_at_0;
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_reduce_degree_success_and_error() -> AlgebraResult<()> {
        // Construct q of degree 3 and elevate once to p of degree 4
        let q = BernsteinPolynomial::new(vec![
            EFloat64::from(1.0),
            EFloat64::from(2.0),
            EFloat64::from(3.0),
            EFloat64::from(4.0),
        ]);
        let p = q.clone().elevate_degree(1);

        // Reducing p by one should exactly recover q
        let reduced = p.reduce_degree(1)?;
        assert_eq!(reduced, q);

        // A random polynomial usually is not exactly reducible
        // Choose coefficients that violate the elevation consistency:
        // For degree 2, reducible iff a2 == 2*a1 - a0. Here: 3 != 2*1 - 0
        let not_square = BernsteinPolynomial::new(vec![
            EFloat64::from(0.0),
            EFloat64::from(1.0),
            EFloat64::from(3.0),
        ]);
        assert!(not_square.reduce_degree(1).is_err());

        Ok(())
    }

    fn ph_cubic_from_w(
        p0: (f64, f64),
        w0: (f64, f64),
        w1: (f64, f64),
    ) -> AlgebraResult<BernsteinPolynomial<Point>> {
        let (a0, b0) = w0;
        let (a1, b1) = w1;

        // complex square and product in R^2
        let sq = |a: f64, b: f64| (a * a - b * b, 2.0 * a * b);
        let mul = |a: f64, b: f64, c: f64, d: f64| (a * c - b * d, a * d + b * c);

        let d0 = sq(a0, b0);
        let d1 = mul(a0, b0, a1, b1);
        let d2 = sq(a1, b1);

        let n = EFloat64::from(3.0); // cubic => divide D_i by 3
        let p0p = Point::from_f64(p0.0, p0.1, 0.0);
        let p1p = p0p + (Point::from_f64(d0.0, d0.1, 0.0) / n)?;
        let p2p = p1p + (Point::from_f64(d1.0, d1.1, 0.0) / n)?;
        let p3p = p2p + (Point::from_f64(d2.0, d2.1, 0.0) / n)?;

        Ok(BernsteinPolynomial::new(vec![p0p, p1p, p2p, p3p]))
    }

    #[test]
    fn test_ph_cubic_derivative_speed_sqrt() -> AlgebraResult<()> {
        // Any non-degenerate choice works; keep z=0
        let c = ph_cubic_from_w((0.0, 0.0), (1.0, 2.0), (2.0, -1.0))?;
        let deriv = c.derivative(); // degree 2
        let speed_sq = deriv.dot(&deriv); // degree 4
        let speed = speed_sq.sqrt().expect("PH curve: sqrt must exist");
        assert_eq!(speed.degree(), 2); // sqrt degree n-1 = 2
        // sanity: (sqrt)^2 == speed_sq
        let back = speed.clone() * speed.clone();
        assert!(
            back.coefficients
                .iter()
                .zip(speed_sq.coefficients.iter())
                .all(|(a, b)| *a == *b)
        );

        // Render the curve and the normalized derivative using the speed polynomial

        // Sample points along the curve
        let n_samples = 32;

        // Compose the scene
        let mut scene = PrimitiveScene::new();
        scene.add_curve(&c, Color10::Blue)?;
        scene.add_curve(&deriv, Color10::Red)?;

        // render the control polygon of deriv
        for c in deriv.coefficients.iter() {
            scene.add_point(*c, Color10::Green);
        }

        for i in 0..=n_samples {
            let t = EFloat64::from(i as f64 / n_samples as f64);
            let pt = c.eval(t);
            let deriv_vec = deriv.eval(t);
            let speed_val = speed.eval(t);

            // Normalize the derivative (tangent)
            let tangent = if speed_val != EFloat64::zero() {
                (deriv_vec / speed_val)?
            } else {
                Point::zero()
            };

            // For visualization, scale the tangent for clarity
            let tangent_tip = pt + tangent * EFloat64::from(0.5);

            let line = Line::try_new(pt, tangent_tip)?;
            scene.add_line(line, Color10::Red);
        }

        scene.save_to_file("test_outputs/ph_cubic_derivative_speed_sqrt.html")?;

        Ok(())
    }

    fn ph_quintic_from_w(
        p0: (f64, f64),
        w0: (f64, f64),
        w1: (f64, f64),
        w2: (f64, f64),
    ) -> AlgebraResult<BernsteinPolynomial<Point>> {
        let (a0, b0) = w0;
        let (a1, b1) = w1;
        let (a2, b2) = w2;

        // complex square and product in R^2
        let sq = |a: f64, b: f64| (a * a - b * b, 2.0 * a * b);
        let mul = |a: f64, b: f64, c: f64, d: f64| (a * c - b * d, a * d + b * c);

        // Hodograph control points D_k for degree 4 derivative as w^2, with Bernstein weights
        // k = 0: D0 = w0^2
        let d0 = sq(a0, b0);
        // k = 1: D1 = w0 * w1
        let d1 = mul(a0, b0, a1, b1);
        // k = 2: D2 = (1/3) w0 w2 + (2/3) w1^2
        let w0w2 = mul(a0, b0, a2, b2);
        let w1sq = sq(a1, b1);
        let d2 = ((w0w2.0 + 2.0 * w1sq.0) / 3.0, (w0w2.1 + 2.0 * w1sq.1) / 3.0);
        // k = 3: D3 = w1 * w2
        let d3 = mul(a1, b1, a2, b2);
        // k = 4: D4 = w2^2
        let d4 = sq(a2, b2);

        let n = EFloat64::from(5.0); // quintic => divide D_i by 5

        let p0p = Point::from_f64(p0.0, p0.1, 0.0);
        let p1p = p0p + (Point::from_f64(d0.0, d0.1, 0.0) / n)?;
        let p2p = p1p + (Point::from_f64(d1.0, d1.1, 0.0) / n)?;
        let p3p = p2p + (Point::from_f64(d2.0, d2.1, 0.0) / n)?;
        let p4p = p3p + (Point::from_f64(d3.0, d3.1, 0.0) / n)?;
        let p5p = p4p + (Point::from_f64(d4.0, d4.1, 0.0) / n)?;

        Ok(BernsteinPolynomial::new(vec![p0p, p1p, p2p, p3p, p4p, p5p]))
    }

    #[test]
    fn test_ph_quintic_derivative_speed_sqrt() -> AlgebraResult<()> {
        // Construct a PH quintic via quadratic complex hodograph w(t)
        let c = ph_quintic_from_w((0.0, 0.0), (1.0, 2.0), (2.0, -1.0), (1.0, 1.0))?;

        let deriv = c.derivative(); // degree 4
        let speed_sq = deriv.dot(&deriv); // degree 8
        let speed = speed_sq.sqrt().expect("PH curve: sqrt must exist");
        assert_eq!(speed.degree(), 4); // sqrt degree n-1 = 4

        // sanity: (sqrt)^2 == speed_sq
        let back = speed.clone() * speed.clone();
        assert!(
            back.coefficients
                .iter()
                .zip(speed_sq.coefficients.iter())
                .all(|(a, b)| *a == *b)
        );

        // Render the curve and the normalized derivative using the speed polynomial
        let n_samples = 32;
        let mut scene = PrimitiveScene::new();
        scene.add_curve(&c, Color10::Blue)?;
        scene.add_curve(&deriv, Color10::Red)?;

        // render the control polygon of deriv
        for cp in deriv.coefficients.iter() {
            scene.add_point(*cp, Color10::Green);
        }

        for i in 0..=n_samples {
            let t = EFloat64::from(i as f64 / n_samples as f64);
            let pt = c.eval(t);
            let deriv_vec = deriv.eval(t);
            let speed_val = speed.eval(t);

            let tangent = if speed_val != EFloat64::zero() {
                (deriv_vec / speed_val)?
            } else {
                Point::zero()
            };

            let tangent_tip = pt + tangent * EFloat64::from(0.5);
            let line = Line::try_new(pt, tangent_tip)?;
            scene.add_line(line, Color10::Red);
        }

        scene.save_to_file("test_outputs/ph_quintic_derivative_speed_sqrt.html")?;

        Ok(())
    }
}
