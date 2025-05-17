use super::{color::Color10, line::Line, point::Point, triangle::TriangleFace};

pub struct PrimitiveScene {
    pub points: Vec<(Point, Color10)>,
    pub lines: Vec<(Line, Color10)>,
    pub triangles: Vec<(TriangleFace, Color10)>,
}
