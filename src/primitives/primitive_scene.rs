use crate::{
    algebra_error::AlgebraResult, curves::curve_like::CurveLike, rasterize::curve::rasterize_curve,
    renderer::render_scene,
};

use super::{color::Color10, line::Line, point::Point, triangle::TriangleFace};

pub struct PrimitiveScene {
    pub points: Vec<(Point, Color10)>,
    pub lines: Vec<(Line, Color10)>,
    pub triangles: Vec<(TriangleFace, Color10)>,
}

impl PrimitiveScene {
    pub fn new() -> Self {
        Self {
            points: vec![],
            lines: vec![],
            triangles: vec![],
        }
    }

    pub fn add_point(&mut self, point: Point, color: Color10) {
        self.points.push((point, color));
    }

    pub fn add_line(&mut self, line: Line, color: Color10) {
        self.lines.push((line, color));
    }

    pub fn add_triangle(&mut self, triangle: TriangleFace, color: Color10) {
        self.triangles.push((triangle, color));
    }

    pub fn add_curve(&mut self, curve: &dyn CurveLike, color: Color10) -> AlgebraResult<()> {
        let lines = rasterize_curve(curve)?;
        self.lines
            .extend(lines.into_iter().map(|line| (line, color.clone())));
        Ok(())
    }

    pub fn save_to_file(&self, filename: &str) -> AlgebraResult<()> {
        render_scene(self, filename)
    }
}
