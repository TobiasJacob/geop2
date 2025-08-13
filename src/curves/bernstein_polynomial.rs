use std::fmt::Display;
use std::ops::{Add, Mul, Sub};

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

pub trait Zero {
    fn zero() -> Self;
}

impl Zero for EFloat64 {
    fn zero() -> Self {
        EFloat64::zero()
    }
}

impl Zero for Point {
    fn zero() -> Self {
        Point::zero()
    }
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

impl<T> BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
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

// Utility function for binomial coefficients
fn binomial_coefficient(n: usize, k: usize) -> usize {
    if k > n {
        0
    } else {
        (1..=k).fold(1, |acc, i| acc * (n + 1 - i) / i)
    }
}

fn equalize_degree<T>(
    lhs: BernsteinPolynomial<T>,
    rhs: BernsteinPolynomial<T>,
) -> (BernsteinPolynomial<T>, BernsteinPolynomial<T>)
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    let n = lhs.degree();
    let m = rhs.degree();
    if n == m {
        return (lhs, rhs);
    }
    if n < m {
        let lhs = lhs.elevate_degree(m - n);
        (lhs, rhs)
    } else {
        let rhs = rhs.elevate_degree(n - m);
        (lhs, rhs)
    }
}

impl<T> Add for BernsteinPolynomial<T>
where
    T: Zero + Clone + Add<Output = T> + Sub<Output = T> + Mul<EFloat64, Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let (lhs, rhs) = equalize_degree(self, rhs);
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
        let (lhs, rhs) = equalize_degree(self, rhs);
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
impl Mul for BernsteinPolynomial<EFloat64> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
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
                acc = acc + self.coefficients[i] * rhs.coefficients[j] * factor;
            }
            coefficients[k] = acc;
        }

        Self::new(coefficients)
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
}
