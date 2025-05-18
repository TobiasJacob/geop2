use crate::primitives::{
    color::Color10, convex_hull::ConvexHull, line::Line, primitive_scene::PrimitiveScene,
};

pub fn rasterize_convex_hull(convex_hull: &ConvexHull, color: Color10) -> PrimitiveScene {
    let mut scene = PrimitiveScene::new();
    match convex_hull {
        ConvexHull::Point(point) => {
            scene.add_point(*point, color);
        }
        ConvexHull::Line(line) => {
            scene.add_line(*line, color);
        }
        ConvexHull::Triangle(triangle) => {
            scene.add_line(Line::try_new(triangle.a, triangle.b).unwrap(), color);
            scene.add_line(Line::try_new(triangle.b, triangle.c).unwrap(), color);
            scene.add_line(Line::try_new(triangle.c, triangle.a).unwrap(), color);
        }
        ConvexHull::Polyhedron(polyhedron) => {
            for triangle in polyhedron {
                scene.add_line(Line::try_new(triangle.a, triangle.b).unwrap(), color);
                scene.add_line(Line::try_new(triangle.b, triangle.c).unwrap(), color);
                scene.add_line(Line::try_new(triangle.c, triangle.a).unwrap(), color);
            }
        }
    }
    scene
}
