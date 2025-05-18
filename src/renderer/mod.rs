use crate::algebra_error::{AlgebraError, AlgebraResult};
use crate::primitives::{color::Color10, primitive_scene::PrimitiveScene};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

impl From<std::io::Error> for AlgebraError {
    fn from(error: std::io::Error) -> Self {
        AlgebraError::new(error.to_string())
    }
}

// Convert Color10 to hex color string
fn color_to_hex(color: &Color10) -> String {
    match color {
        Color10::Blue => "#1f77b4".to_string(),
        Color10::Orange => "#ff7f0e".to_string(),
        Color10::Green => "#2ca02c".to_string(),
        Color10::Red => "#d62728".to_string(),
        Color10::Purple => "#9467bd".to_string(),
        Color10::Brown => "#8c564b".to_string(),
        Color10::Pink => "#e377c2".to_string(),
        Color10::Gray => "#7f7f7f".to_string(),
        Color10::Olive => "#bcbd22".to_string(),
        Color10::Cyan => "#17becf".to_string(),
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

    // Convert points to JavaScript array with colors
    let points_data: Vec<(f64, f64, f64, String)> = scene
        .points
        .iter()
        .map(|(p, color)| {
            (
                p.x.lower_bound,
                p.y.lower_bound,
                p.z.lower_bound,
                color_to_hex(color),
            )
        })
        .collect();

    // Convert lines to JavaScript array with colors
    let lines_data: Vec<(f64, f64, f64, f64, f64, f64, String)> = scene
        .lines
        .iter()
        .map(|(l, color)| {
            let start = l.start();
            let end = l.end();
            (
                start.x.to_f64(),
                start.y.to_f64(),
                start.z.to_f64(),
                end.x.to_f64(),
                end.y.to_f64(),
                end.z.to_f64(),
                color_to_hex(color),
            )
        })
        .collect();

    // Convert triangles to JavaScript array with colors
    let triangles_data: Vec<(f64, f64, f64, f64, f64, f64, f64, f64, f64, String)> = scene
        .triangles
        .iter()
        .map(|(t, color)| {
            (
                t.a.x.to_f64(),
                t.a.y.to_f64(),
                t.a.z.to_f64(),
                t.b.x.to_f64(),
                t.b.y.to_f64(),
                t.b.z.to_f64(),
                t.c.x.to_f64(),
                t.c.y.to_f64(),
                t.c.z.to_f64(),
                color_to_hex(color),
            )
        })
        .collect();

    // Replace the placeholders with actual data
    let template = template.replace(
        "/* Add points here */",
        &format!(
            "{}",
            points_data
                .iter()
                .map(|(x, y, z, color)| format!("[{},{},{},'{}']", x, y, z, color))
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
                .map(|(x1, y1, z1, x2, y2, z2, color)| {
                    format!("[{},{},{},{},{},{},'{}']", x1, y1, z1, x2, y2, z2, color)
                })
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
                .map(|(x1, y1, z1, x2, y2, z2, x3, y3, z3, color)| {
                    format!(
                        "[{},{},{},{},{},{},{},{},{},'{}']",
                        x1, y1, z1, x2, y2, z2, x3, y3, z3, color
                    )
                })
                .collect::<Vec<_>>()
                .join(",")
        ),
    );

    // Replace the debug text placeholder
    let template = template.replace("Debug text goes here...", &scene.debug_text);

    // Write the modified template to the file
    writeln!(file, "{}", template)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::primitives::{color::Color10, line::Line, point::Point, triangle::TriangleFace};

    use super::*;
    use std::fs;

    #[test]
    fn test_render_scene() -> AlgebraResult<()> {
        // Create a test scene with various primitives
        let scene = PrimitiveScene {
            // Add points in a cube pattern
            points: vec![
                (Point::from_f64(0.0, 0.0, 0.0), Color10::Red), // origin
                (Point::from_f64(1.0, 0.0, 0.0), Color10::Gray), // x-axis
                (Point::from_f64(0.0, 1.0, 0.0), Color10::Gray), // y-axis
                (Point::from_f64(0.0, 0.0, 1.0), Color10::Gray), // z-axis
            ],
            // Add lines forming a cube
            lines: vec![
                // Bottom face
                (
                    Line::try_new(
                        Point::from_f64(0.0, 0.0, 0.0),
                        Point::from_f64(1.0, 0.0, 0.0),
                    )
                    .unwrap(),
                    Color10::Blue,
                ),
                (
                    Line::try_new(
                        Point::from_f64(1.0, 0.0, 0.0),
                        Point::from_f64(1.0, 1.0, 0.0),
                    )
                    .unwrap(),
                    Color10::Orange,
                ),
                (
                    Line::try_new(
                        Point::from_f64(1.0, 1.0, 0.0),
                        Point::from_f64(0.0, 1.0, 0.0),
                    )
                    .unwrap(),
                    Color10::Green,
                ),
                (
                    Line::try_new(
                        Point::from_f64(0.0, 1.0, 0.0),
                        Point::from_f64(0.0, 0.0, 0.0),
                    )
                    .unwrap(),
                    Color10::Red,
                ),
                // Top face
                (
                    Line::try_new(
                        Point::from_f64(0.0, 0.0, 1.0),
                        Point::from_f64(1.0, 0.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Purple,
                ),
                (
                    Line::try_new(
                        Point::from_f64(1.0, 0.0, 1.0),
                        Point::from_f64(1.0, 1.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Brown,
                ),
                (
                    Line::try_new(
                        Point::from_f64(1.0, 1.0, 1.0),
                        Point::from_f64(0.0, 1.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Pink,
                ),
                (
                    Line::try_new(
                        Point::from_f64(0.0, 1.0, 1.0),
                        Point::from_f64(0.0, 0.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Blue,
                ),
                // Connecting lines
                (
                    Line::try_new(
                        Point::from_f64(0.0, 0.0, 0.0),
                        Point::from_f64(0.0, 0.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Orange,
                ),
                (
                    Line::try_new(
                        Point::from_f64(1.0, 0.0, 0.0),
                        Point::from_f64(1.0, 0.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Green,
                ),
                (
                    Line::try_new(
                        Point::from_f64(1.0, 1.0, 0.0),
                        Point::from_f64(1.0, 1.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Red,
                ),
                (
                    Line::try_new(
                        Point::from_f64(0.0, 1.0, 0.0),
                        Point::from_f64(0.0, 1.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Purple,
                ),
            ],
            // Add triangles forming a tetrahedron
            triangles: vec![
                // Base triangle
                (
                    TriangleFace::try_new(
                        Point::from_f64(0.0, 0.0, 0.0),
                        Point::from_f64(1.0, 0.0, 0.0),
                        Point::from_f64(0.0, 1.0, 0.0),
                    )
                    .unwrap(),
                    Color10::Blue,
                ),
                // Side triangles
                (
                    TriangleFace::try_new(
                        Point::from_f64(0.0, 0.0, 0.0),
                        Point::from_f64(0.0, 1.0, 0.0),
                        Point::from_f64(0.0, 0.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Orange,
                ),
                (
                    TriangleFace::try_new(
                        Point::from_f64(0.0, 0.0, 0.0),
                        Point::from_f64(0.0, 0.0, 1.0),
                        Point::from_f64(1.0, 0.0, 0.0),
                    )
                    .unwrap(),
                    Color10::Green,
                ),
                (
                    TriangleFace::try_new(
                        Point::from_f64(0.0, 1.0, 0.0),
                        Point::from_f64(1.0, 0.0, 0.0),
                        Point::from_f64(0.0, 0.0, 1.0),
                    )
                    .unwrap(),
                    Color10::Red,
                ),
            ],
            debug_text: "Test debug text".to_string(),
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
