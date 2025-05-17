use super::{line::Line, point::Point, triangle::TriangleFace};

pub struct PrimitiveScene {
    pub points: Vec<Point>,
    pub lines: Vec<Line>,
    pub triangles: Vec<TriangleFace>,
}
