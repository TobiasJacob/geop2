pub struct Face {
    pub surface: Box<dyn SurfaceLike>,
    pub curves: Vec<Box<dyn CurveLike>>,
}

impl Face {
    pub fn get_convex_hull(&self) -> GeometryResult<BoundingBox> {
        let mut convex_hull = self.surface.get_convex_hull()?;
        for curve in self.curves.iter() {
            convex_hull.merge(&curve.get_convex_hull()?);
        }
        Ok(convex_hull)
    }
}
