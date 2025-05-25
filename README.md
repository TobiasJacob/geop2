# Geop 2

Geop 2 is an open source CAD kernel. It is the successor to [Geop](https://github.com/TobiasJacob/geop), which worked as a proof of concept. Soon, we will also add features from [Isotope](https://github.com/CADmium-Co/ISOtope), a high-performance sketch constraint solver. We use some lessons learned from Geop to build a more robust and feature-rich kernel. Here are some of the key lessons learned:
- We ditched math heavy algorithms in favor of subdivions schemes.
- We are NURBS first. Once intersections, boolean operations, fillets, etc. are working well for NURBS, we can convert other geometric entities like planes or spheres to nurbs.
- For debugging, we are still rendering everything into a book, but now we render to html instead of pngs.

## Roadmap
- [x] NURBS curves and surfaces
- [x] Convex hulls decomposition and intersection
- [x] Subdivion schemes
    - [x] Curve subdivision
    - [x] Surface subdivision
- [ ] Intersections
    - [x] Curve-curve intersections
    - [ ] Curve-surface intersections
    - [ ] Surface-surface intersections
- [ ] Boolean operations
- [ ] Fillets and Chamfers
- [ ] Sketch constraint solver
- [ ] Export to STEP
- [ ] Export to STL

- NURBS curves and surfaces
- Boolean operations

## Demos

![Convex Hull](docs\ConvexHull.png)

![Curve Intersection](docs\Intersection.png)
