use crate::{
    algebra_error::AlgebraResult,
    curves::curve_like::CurveLike,
    face::Face,
    rasterize::{
        convex_hull::rasterize_convex_hull, curve::rasterize_curve, surface::rasterize_surface,
    },
    renderer::render_scene,
};

use super::{
    color::Color10, convex_hull::ConvexHull, line::Line, point::Point, triangle::TriangleFace,
};

#[derive(Debug, Clone)]
pub struct PrimitiveScene {
    pub points: Vec<(Point, Color10)>,
    pub lines: Vec<(Line, Color10)>,
    pub triangles: Vec<(TriangleFace, Color10)>,
    pub debug_text: String,
}

impl PrimitiveScene {
    pub fn new() -> Self {
        Self {
            points: vec![],
            lines: vec![],
            triangles: vec![],
            debug_text: String::new(),
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

    pub fn set_debug_text(&mut self, text: String) {
        self.debug_text = text;
    }

    pub fn add_debug_text(&mut self, text: String) {
        if !self.debug_text.is_empty() {
            self.debug_text.push('\n');
        }
        self.debug_text.push_str(&text);
    }

    pub fn add_scene(&mut self, scene: PrimitiveScene) {
        self.points.extend(scene.points);
        self.lines.extend(scene.lines);
        self.triangles.extend(scene.triangles);
        if !scene.debug_text.is_empty() {
            if !self.debug_text.is_empty() {
                self.debug_text.push('\n');
            }
            self.debug_text.push_str(&scene.debug_text);
        }
    }

    pub fn add_convex_hull(&mut self, convex_hull: ConvexHull, color: Color10) {
        let scene = rasterize_convex_hull(&convex_hull, color);
        self.add_scene(scene);
    }

    pub fn add_curve(&mut self, curve: &dyn CurveLike, color: Color10) -> AlgebraResult<()> {
        let lines = rasterize_curve(curve)?;
        self.lines
            .extend(lines.into_iter().map(|line| (line, color.clone())));
        Ok(())
    }

    pub fn add_face(&mut self, face: &Face, color: Color10) -> AlgebraResult<()> {
        let triangles = rasterize_surface(face)?;
        self.triangles.extend(
            triangles
                .into_iter()
                .map(|triangle| (triangle, color.clone())),
        );
        Ok(())
    }

    pub fn save_to_file(&self, filename: &str) -> AlgebraResult<()> {
        render_scene(self, filename)
    }
}
