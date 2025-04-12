// Nurbs implement functions like transform, revolve, extrude, etc.

// CurveLike is a trait that all curves should implement.
// For circles we use a trick. A circle always has an anchor point, such that we have a defined start and end.
pub trait CurveLike {
    // Transform
    fn transform(&self, transform: Transform) -> Self;

    // Nurb
    fn to_nurbs(&self) -> Nurbs;

    // Change the direction of the curve
    fn neg(&self) -> Curve;

    // Normalized Tangent / Direction of the curve at the given point.
    fn tangent(&self, p: Point) -> GeometryResult<Point>;

    // Curvature of the curve at the given point.
    fn curvature(&self, p: Point) -> GeometryResult<EFloat64>;

    // Checks if point is on the curve. (between start and end)
    fn on_curve(&self, p: Point) -> bool;

    // Returns the distance between x and y. Fails if x and y are not on the curve.
    fn distance(&self, x: Point, y: Point) -> GeometryResult<EFloat64>;

    // Interpolate between start and end at t. t is between 0 and 1. Keep in mind that it might be a line or a circle.
    // This is the only function exposing a curve parameter to the outside world.
    // fn interpolate(&self, t: f64) -> GeometryResult<Point>;

    // Rasterize the curve into n points.
    fn rasterize(&self, n: usize) -> GeometryResult<Vec<Point>>;

    // Finds the closest point on the curve to the given point.
    fn project(&self, p: Point) -> Point;

    // Get the midpoint between start and end.
    // This will guarantee that between(start, midpoint, end) is true and midpoint != start and midpoint != end.
    // If start or end is None, the midpoint is a point that is a unit distance away from the other point.
    // This operation will fail if start == end.
    fn get_inner_point(&self) -> GeometryResult<Point>;

    // Split at point p.
    fn split(&self, p: Point) -> (Self, Self);

    // Subdivide the curve at t.
    fn subdivide(&self) -> (Self, Self);

    // Returns a bounding box that contains the curve.
    fn get_convex_hull(&self) -> GeometryResult<BoundingBox>;
}
