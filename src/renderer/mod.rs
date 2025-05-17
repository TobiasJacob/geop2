use crate::algebra_error::{AlgebraError, AlgebraResult};
use crate::primitives::primitive_scene::PrimitiveScene;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

impl From<std::io::Error> for AlgebraError {
    fn from(error: std::io::Error) -> Self {
        AlgebraError::new(error.to_string())
    }
}

// Render a scene into an html file
pub fn render_scene(scene: &PrimitiveScene, file_path: &str) -> AlgebraResult<()> {
    // if the directory does not exist, create it
    let dir = Path::new(file_path).parent().unwrap();
    if !dir.exists() {
        fs::create_dir_all(dir)?;
    }
    let mut file = File::create(file_path)?;

    // Read the template
    let template = include_str!("template.html");

    // Convert points to JavaScript array
    let points_data: Vec<f64> = scene
        .points
        .iter()
        .flat_map(|p| vec![p.x.lower_bound, p.y.lower_bound, p.z.lower_bound])
        .collect();

    // Convert lines to JavaScript array
    let lines_data: Vec<f64> = scene
        .lines
        .iter()
        .flat_map(|l| {
            let start = l.start();
            let end = l.end();
            vec![
                start.x.to_f64(),
                start.y.to_f64(),
                start.z.to_f64(),
                end.x.to_f64(),
                end.y.to_f64(),
                end.z.to_f64(),
            ]
        })
        .collect();

    // Convert triangles to JavaScript array
    let triangles_data: Vec<f64> = scene
        .triangles
        .iter()
        .flat_map(|t| {
            vec![
                t.a.x.to_f64(),
                t.a.y.to_f64(),
                t.a.z.to_f64(),
                t.b.x.to_f64(),
                t.b.y.to_f64(),
                t.b.z.to_f64(),
                t.c.x.to_f64(),
                t.c.y.to_f64(),
                t.c.z.to_f64(),
            ]
        })
        .collect();

    // Replace the placeholders with actual data
    let template = template.replace(
        "/* Add points here */",
        &format!(
            "{}",
            points_data
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let template = template.replace(
        "/* Add lines here */",
        &format!(
            "{}",
            lines_data
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let template = template.replace(
        "/* Add triangles here */",
        &format!(
            "{}",
            triangles_data
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
    );

    // Write the modified template to the file
    writeln!(file, "{}", template)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::primitives::{line::Line, point::Point, triangle::TriangleFace};

    use super::*;
    use std::fs;

    #[test]
    fn test_render_scene() -> AlgebraResult<()> {
        // Create a test scene with various primitives
        let scene = PrimitiveScene {
            // Add points in a cube pattern
            points: vec![
                Point::from_f64(0.0, 0.0, 0.0), // origin
                Point::from_f64(1.0, 0.0, 0.0), // x-axis
                Point::from_f64(0.0, 1.0, 0.0), // y-axis
                Point::from_f64(0.0, 0.0, 1.0), // z-axis
            ],
            // Add lines forming a cube
            lines: vec![
                // Bottom face
                Line::try_new(
                    Point::from_f64(0.0, 0.0, 0.0),
                    Point::from_f64(1.0, 0.0, 0.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(1.0, 0.0, 0.0),
                    Point::from_f64(1.0, 1.0, 0.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(1.0, 1.0, 0.0),
                    Point::from_f64(0.0, 1.0, 0.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(0.0, 1.0, 0.0),
                    Point::from_f64(0.0, 0.0, 0.0),
                )
                .unwrap(),
                // Top face
                Line::try_new(
                    Point::from_f64(0.0, 0.0, 1.0),
                    Point::from_f64(1.0, 0.0, 1.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(1.0, 0.0, 1.0),
                    Point::from_f64(1.0, 1.0, 1.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(1.0, 1.0, 1.0),
                    Point::from_f64(0.0, 1.0, 1.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(0.0, 1.0, 1.0),
                    Point::from_f64(0.0, 0.0, 1.0),
                )
                .unwrap(),
                // Connecting lines
                Line::try_new(
                    Point::from_f64(0.0, 0.0, 0.0),
                    Point::from_f64(0.0, 0.0, 1.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(1.0, 0.0, 0.0),
                    Point::from_f64(1.0, 0.0, 1.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(1.0, 1.0, 0.0),
                    Point::from_f64(1.0, 1.0, 1.0),
                )
                .unwrap(),
                Line::try_new(
                    Point::from_f64(0.0, 1.0, 0.0),
                    Point::from_f64(0.0, 1.0, 1.0),
                )
                .unwrap(),
            ],
            // Add triangles forming a tetrahedron
            triangles: vec![
                // Base triangle
                TriangleFace::try_new(
                    Point::from_f64(0.0, 0.0, 0.0),
                    Point::from_f64(1.0, 0.0, 0.0),
                    Point::from_f64(0.0, 1.0, 0.0),
                )
                .unwrap(),
                // Side triangles
                TriangleFace::try_new(
                    Point::from_f64(0.0, 0.0, 0.0),
                    Point::from_f64(0.0, 1.0, 0.0),
                    Point::from_f64(0.0, 0.0, 1.0),
                )
                .unwrap(),
                TriangleFace::try_new(
                    Point::from_f64(0.0, 0.0, 0.0),
                    Point::from_f64(0.0, 0.0, 1.0),
                    Point::from_f64(1.0, 0.0, 0.0),
                )
                .unwrap(),
                TriangleFace::try_new(
                    Point::from_f64(0.0, 1.0, 0.0),
                    Point::from_f64(1.0, 0.0, 0.0),
                    Point::from_f64(0.0, 0.0, 1.0),
                )
                .unwrap(),
            ],
        };

        // Render the scene to a test file
        let test_file = "test_outputs/test_scene.html";
        render_scene(&scene, test_file)?;

        // Verify the file was created and contains the expected content
        let content = fs::read_to_string(test_file)?;
        assert!(content.contains("THREE.Scene"));
        assert!(content.contains("THREE.Points"));
        assert!(content.contains("THREE.LineSegments"));
        assert!(content.contains("THREE.Mesh"));

        // Clean up the test file
        // fs::remove_file(test_file)?;

        Ok(())
    }
}
