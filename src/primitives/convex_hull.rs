use std::fmt::Display;

use crate::algebra_error::{AlgebraError, AlgebraResult, WithContext};

use crate::primitives::intersecting::{
    line::line_line_intersection,
    line::line_triangle_intersection,
    point::point_line_intersection,
    point::point_triangle_intersection,
    polyhedron::{
        line_polyhedron_intersection, point_polyhedron_intersection,
        polyhedron_polyhedron_intersection, triangle_polyhedron_intersection,
    },
    triangle::triangle_triangle_intersection,
};
use crate::primitives::line::Line;
use crate::primitives::point::Point;
use crate::primitives::triangle::TriangleFace;

use super::efloat::EFloat64;

/// Computes the distance from a point `p` to the line defined by points `a` and `b`.
pub fn distance_point_to_line(p: &Point, a: &Point, b: &Point) -> EFloat64 {
    let ap = *p - *a;
    let ab = *b - *a;
    let ab_norm = ab.norm();
    if ab_norm == 0.0 {
        return ap.norm();
    }
    let cross = ap.cross(ab);
    let cross_norm = cross.norm();
    (cross_norm / ab_norm).unwrap()
}

/// A convex hull that can represent different geometric objects based on the number of points:
/// - Single point
/// - Line segment (two points)
/// - Triangle (three points)
/// - Convex polyhedron (more than three points)
#[derive(Debug, Clone)]
pub enum ConvexHull {
    Point(Point),
    Line(Line),
    Triangle(TriangleFace),
    Polyhedron(Vec<TriangleFace>),
}

impl ConvexHull {
    /// Creates a new convex hull from a set of points using the Quickhull algorithm.
    pub fn try_new(points: Vec<Point>) -> AlgebraResult<Self> {
        let context = |err: AlgebraError| {
            let mut message = format!("Creating convex hull from points");
            // now print each point to the message
            for point in points.iter() {
                message.push_str(&format!(" {}", point));
            }
            err.with_context(message)
        };

        // Filter out duplicate points
        let mut unique_points = Vec::new();
        for point in points.clone() {
            if !unique_points.iter().any(|&p| p == point) {
                unique_points.push(point);
            }
        }

        match unique_points.len() {
            0 => Err(AlgebraError::new(
                "Cannot create convex hull from empty set of points".to_string(),
            )),
            1 => Ok(Self::Point(unique_points[0])),
            2 => Ok(Self::Line(Line::try_new(
                unique_points[0],
                unique_points[1],
            )?)),
            3 => {
                // Make sure the points are not collinear.
                let cross = (unique_points[1] - unique_points[0])
                    .cross(unique_points[2] - unique_points[0]);
                if cross.norm() == 0.0 {
                    // check which one is the farthest from the other two
                    let d1 = (unique_points[0] - unique_points[1]).norm();
                    let d2 = (unique_points[0] - unique_points[2]).norm();
                    let d3 = (unique_points[1] - unique_points[2]).norm();
                    if d1 > d2 && d1 > d3 {
                        unique_points = vec![unique_points[0], unique_points[2]];
                    } else if d2 > d3 {
                        unique_points = vec![unique_points[0], unique_points[1]];
                    } else {
                        unique_points = vec![unique_points[1], unique_points[2]];
                    }
                    return Ok(Self::Line(Line::try_new(
                        unique_points[0],
                        unique_points[1],
                    )?));
                } else {
                    Ok(Self::Triangle(
                        TriangleFace::try_new(unique_points[0], unique_points[1], unique_points[2])
                            .with_context(&context)?,
                    ))
                }
            }
            _ => {
                // 1. Determine two unique_points with the minimum and maximum x-values (they are definitely part of the hull).
                let mut min_x = unique_points[0];
                let mut max_x = unique_points[0];
                for p in unique_points.iter() {
                    if p.x < min_x.x {
                        min_x = *p;
                    }
                    if p.x > max_x.x {
                        max_x = *p;
                    }
                }

                // 2. Determine the point that is farthest from the line (min_x, max_x).
                let mut max_distance = EFloat64::from(-1.0);
                let mut third_point = unique_points[0];
                for p in unique_points.iter() {
                    let d = distance_point_to_line(p, &min_x, &max_x);
                    if d > max_distance {
                        max_distance = d;
                        third_point = *p;
                    }
                }

                // 3. Determine the point that is farthest from the plane (min_x, max_x, third_point).
                let temp_face = TriangleFace::try_new(min_x, max_x, third_point)?;
                max_distance = EFloat64::from(0.0);

                // 3.1 Start with an initial point that is different from min_x, max_x, and third_point. This is imporant in the coplanar case.
                let mut fourth_point = unique_points[0];
                for p in unique_points.iter() {
                    if p != &min_x && p != &max_x && p != &third_point {
                        fourth_point = *p;
                        break;
                    }
                }

                for p in unique_points.iter() {
                    let d = temp_face.distance_to_point(p).abs();
                    if d > max_distance {
                        max_distance = d;
                        fourth_point = *p;
                    }
                }

                // 4. Create the initial tetrahedron (four faces).
                let mut faces = Vec::new();
                // Check the orientation of the fourth point relative to the plane (min_x, max_x, third_point).
                if temp_face.distance_to_point(&fourth_point) < 0.0 {
                    faces.push(
                        TriangleFace::try_new(min_x, max_x, third_point).with_context(&context)?,
                    );
                    faces.push(
                        TriangleFace::try_new(min_x, third_point, fourth_point)
                            .with_context(&context)?,
                    );
                    faces.push(
                        TriangleFace::try_new(min_x, fourth_point, max_x).with_context(&context)?,
                    );
                    faces.push(
                        TriangleFace::try_new(max_x, fourth_point, third_point)
                            .with_context(&context)?,
                    );
                } else {
                    faces.push(
                        TriangleFace::try_new(min_x, third_point, max_x).with_context(&context)?,
                    );
                    faces.push(
                        TriangleFace::try_new(min_x, fourth_point, third_point)
                            .with_context(&context)?,
                    );
                    faces.push(
                        TriangleFace::try_new(min_x, max_x, fourth_point).with_context(&context)?,
                    );
                    faces.push(
                        TriangleFace::try_new(max_x, third_point, fourth_point)
                            .with_context(&context)?,
                    );
                }

                // 5. Iterative expansion of the hull:
                // For every point that lies outside (i.e., in front of a face), all visible faces are removed
                // and replaced by new faces connecting the "horizon" (boundary) of the visible faces with the point.
                let mut changed = true;
                while changed {
                    changed = false;
                    while let Some(p) = unique_points.pop() {
                        // Check if a point lies outside of any face.
                        let mut is_outside = false;
                        for face in &faces {
                            if face.distance_to_point(&p) > 0.0 {
                                is_outside = true;
                                break;
                            }
                        }
                        if is_outside {
                            changed = true;
                            // Find all faces that can "see" the point (visible faces).
                            let mut visible_faces = Vec::new();
                            for (i, face) in faces.iter().enumerate() {
                                if face.distance_to_point(&p) > 0.0 {
                                    visible_faces.push(i);
                                }
                            }
                            // Determine the boundary edges (edges that belong only to one visible face).
                            let mut boundary_edges = Vec::new();
                            for &i in &visible_faces {
                                let face = &faces[i];
                                let edges =
                                    vec![(face.a, face.b), (face.b, face.c), (face.c, face.a)];
                                for edge in edges {
                                    // Check if this edge also occurs in a non-visible face.
                                    let mut shared = false;
                                    for (j, other_face) in faces.iter().enumerate() {
                                        if !visible_faces.contains(&j) {
                                            let other_edges = vec![
                                                (other_face.a, other_face.b),
                                                (other_face.b, other_face.c),
                                                (other_face.c, other_face.a),
                                            ];
                                            if other_edges.contains(&(edge.1, edge.0)) {
                                                shared = true;
                                                break;
                                            }
                                        }
                                    }
                                    if shared {
                                        boundary_edges.push(edge);
                                    }
                                }
                            }

                            // Remove all visible faces.
                            faces = faces
                                .into_iter()
                                .enumerate()
                                .filter(|(i, _)| !visible_faces.contains(i))
                                .map(|(_, f)| f)
                                .collect();
                            // Create new faces connecting the boundary edges with the point.
                            for edge in boundary_edges {
                                faces.push(TriangleFace::try_new(edge.0, edge.1, p)?);
                            }
                            // Restart the outer loop after the change.
                            break;
                        }
                    }
                }

                Ok(Self::Polyhedron(faces))
            }
        }
    }

    /// Checks if this convex hull intersects with another convex hull using the separating axis theorem.
    pub fn intersects(&self, other: &ConvexHull) -> bool {
        match (self, other) {
            (Self::Point(p1), Self::Point(p2)) => p1 == p2,
            (Self::Point(p), Self::Line(l)) | (Self::Line(l), Self::Point(p)) => {
                point_line_intersection(p, l)
            }
            (Self::Point(p), Self::Triangle(t)) | (Self::Triangle(t), Self::Point(p)) => {
                point_triangle_intersection(p, t)
            }
            (Self::Point(p), Self::Polyhedron(faces))
            | (Self::Polyhedron(faces), Self::Point(p)) => point_polyhedron_intersection(p, faces),
            (Self::Line(l1), Self::Line(l2)) => line_line_intersection(l1, l2),
            (Self::Line(l), Self::Triangle(t)) | (Self::Triangle(t), Self::Line(l)) => {
                line_triangle_intersection(l, t)
            }
            (Self::Line(l), Self::Polyhedron(faces)) | (Self::Polyhedron(faces), Self::Line(l)) => {
                line_polyhedron_intersection(l, faces)
            }
            (Self::Triangle(t1), Self::Triangle(t2)) => triangle_triangle_intersection(t1, t2),
            (Self::Triangle(t), Self::Polyhedron(faces))
            | (Self::Polyhedron(faces), Self::Triangle(t)) => {
                triangle_polyhedron_intersection(t, faces)
            }
            (Self::Polyhedron(faces1), Self::Polyhedron(faces2)) => {
                polyhedron_polyhedron_intersection(faces1, faces2)
            }
        }
    }

    /// Checks if a point is contained within the convex hull.
    pub fn contains_point(&self, point: &Point) -> bool {
        match self {
            Self::Point(p) => p == point,
            Self::Line(l) => point_line_intersection(point, l),
            Self::Triangle(t) => point_triangle_intersection(point, t),
            Self::Polyhedron(faces) => point_polyhedron_intersection(point, faces),
        }
    }

    pub fn to_points(&self) -> Vec<Point> {
        match self {
            Self::Point(p) => vec![p.clone()],
            Self::Line(l) => vec![l.start().clone(), l.end().clone()],
            Self::Triangle(t) => vec![t.a.clone(), t.b.clone(), t.c.clone()],
            Self::Polyhedron(faces) => faces
                .iter()
                .flat_map(|f| vec![f.a.clone(), f.b.clone(), f.c.clone()])
                .collect(),
        }
    }

    pub fn to_point(&self) -> Point {
        let points = self.to_points();
        let mut result = points[0].clone();
        for p in points {
            result.x = EFloat64::new_union_efloat(result.x, p.x);
            result.y = EFloat64::new_union_efloat(result.y, p.y);
            result.z = EFloat64::new_union_efloat(result.z, p.z);
        }
        result
    }
}

impl Display for ConvexHull {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Point(p) => write!(f, "{}", p),
            Self::Line(l) => write!(f, "{}", l),
            Self::Triangle(t) => write!(f, "{}", t),
            Self::Polyhedron(faces) => {
                for face in faces {
                    writeln!(f, "{}", face)?;
                }
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convex_hull_creation() -> AlgebraResult<()> {
        // Test single point
        let point = Point::from_f64(1.0, 2.0, 3.0);
        let hull = ConvexHull::try_new(vec![point])?;
        assert!(matches!(hull, ConvexHull::Point(p) if p == point));

        // Test line segment
        let p1 = Point::from_f64(0.0, 0.0, 0.0);
        let p2 = Point::from_f64(1.0, 1.0, 1.0);
        let hull = ConvexHull::try_new(vec![p1, p2])?;
        assert!(matches!(hull, ConvexHull::Line(_)));

        // Test triangle
        let p3 = Point::from_f64(0.0, 1.0, 0.0);
        let hull = ConvexHull::try_new(vec![p1, p2, p3])?;
        assert!(matches!(hull, ConvexHull::Triangle(_)));

        // Test polyhedron
        let p4 = Point::from_f64(0.0, 0.0, 1.0);
        let hull = ConvexHull::try_new(vec![p1, p2, p3, p4])?;
        assert!(matches!(hull, ConvexHull::Polyhedron(_)));

        Ok(())
    }

    #[test]
    fn test_convex_hull_intersection() -> AlgebraResult<()> {
        // Create two non-intersecting convex hulls
        let points1 = vec![
            Point::from_f64(0.0, 0.0, 0.0),
            Point::from_f64(1.0, 0.0, 0.0),
            Point::from_f64(0.0, 1.0, 0.0),
            Point::from_f64(0.0, 0.0, 1.0),
        ];
        let points2 = vec![
            Point::from_f64(2.0, 0.0, 0.0),
            Point::from_f64(3.0, 0.0, 0.0),
            Point::from_f64(2.0, 1.0, 0.0),
            Point::from_f64(2.0, 0.0, 1.0),
        ];
        let hull1 = ConvexHull::try_new(points1)?;
        let hull2 = ConvexHull::try_new(points2)?;

        println!("hull1: {}", hull1);
        println!("hull2: {}", hull2);
        assert!(!hull1.intersects(&hull2));

        // Create two intersecting convex hulls
        let points3 = vec![
            Point::from_f64(0.0, 0.0, 0.0),
            Point::from_f64(1.0, 0.0, 0.0),
            Point::from_f64(0.0, 1.0, 0.0),
            Point::from_f64(0.0, 0.0, 1.0),
        ];
        let points4 = vec![
            Point::from_f64(0.5, 0.0, 0.0),
            Point::from_f64(1.5, 0.0, 0.0),
            Point::from_f64(0.5, 1.0, 0.0),
            Point::from_f64(0.5, 0.0, 1.0),
        ];
        let hull3 = ConvexHull::try_new(points3)?;
        let hull4 = ConvexHull::try_new(points4)?;
        assert!(hull3.intersects(&hull4));

        Ok(())
    }
}
