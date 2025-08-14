pub mod algebra_error;
pub mod bernstein;
pub mod contour;
pub mod curves;
pub mod face;
pub mod intersection;
pub mod primitives;
pub mod rasterize;
pub mod renderer;
pub mod surfaces;
pub mod zero;

// Utility function for binomial coefficients
pub fn binomial_coefficient(n: usize, k: usize) -> usize {
    if k > n {
        0
    } else {
        (1..=k).fold(1, |acc, i| acc * (n + 1 - i) / i)
    }
}
