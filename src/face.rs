pub struct Face {
    pub surface: Box<dyn SurfaceLike>,
    pub curves: Vec<Box<dyn CurveLike>>,
}

impl Face {}
