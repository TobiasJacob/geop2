use std::fmt::Display;

use crate::{
    algebra_error::AlgebraResult,
    curves::curve_like::CurveLike,
    primitives::{convex_hull::ConvexHull, efloat::EFloat64, point::Point},
};

// Represents a polynomial in the form of a_{0} B_{0,n}
#[derive(Debug, Clone)]
pub struct BernsteinPolynomial {
    pub coefficients: Vec<Point>,
}

impl BernsteinPolynomial {
    pub fn bernstein_basis(i: usize, n: usize, up_vector: Point, direction_vector: Point) -> Self {
        let mut coefficients = vec![Point::zero(); n + 1];
        for j in 0..=n {
            coefficients[j] = direction_vector.clone() * EFloat64::from(j as f64 / n as f64);
        }
        coefficients[i] = coefficients[i] + up_vector;
        Self::new(coefficients)
    }
}

impl BernsteinPolynomial {
    pub fn new(coefficients: Vec<Point>) -> Self {
        Self { coefficients }
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    //$$ c_i^{n+r} = \sum_{j = max(0, i - r)}^{min(n, i)} \frac{\binom{r}{i - j} \binom{n}{j}}{\binom{n + r}{i}} c_i^n $$
    pub fn elevate_degree(&self, r: usize) -> Self {
        let n = self.degree();
        let mut new_coeffs = vec![Point::zero(); n + r + 1];

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
    pub fn subdivide(&self, t: EFloat64) -> (BernsteinPolynomial, BernsteinPolynomial) {
        let mut beta = self.coefficients.clone();
        let n = beta.len();
        let mut left = vec![Point::zero(); n];
        let mut right = vec![Point::zero(); n];

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

    pub fn eval(&self, t: EFloat64) -> Point {
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

    pub fn derivative(&self) -> BernsteinPolynomial {
        let n = self.degree();
        if n == 0 {
            return BernsteinPolynomial::new(vec![Point::zero()]);
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

impl CurveLike for BernsteinPolynomial {
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

impl Display for BernsteinPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        let n = self.degree();
        for (i, coeff) in self.coefficients.iter().enumerate() {
            if *coeff != Point::zero() {
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
    fn test_bernstein_polynomial() -> BernsteinPolynomial {
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
}
