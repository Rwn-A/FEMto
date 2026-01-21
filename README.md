# FEMto
Finite-Element Method simulation written completely from scratch in Odin.

> [!NOTE]
> Project is in infancy, docs limited for now.

## Goals
FEMto aims to build a lightweight, method-agnostic finite-element core with minimal abstractions.
The core is intended to stand alone for prototyping, while also being the base layer for a higher-level physics interface.
Eventually, we plan to expose the physics interface, and physics built on top of it, through a runtime configuration layer.

## Core Functionality
The core library supports the following:

- Scalar/vector lagrange basis up to P2
- 1st and 2nd order 1D-3D meshes
- Global continuity helpers for the continuous-galerkin method
- Gmsh input and VTU output

the core library consists of the `src/fem` and `src/la` directories, these can be used independantly of the rest of the project.

Performance optimizations have yet to be implemented.

### Mesh Restrictions
The gmsh parser expects the following:

- **All elements must have a physical group or they are ignored**
- Format in **version 2 binary**
- Well-behaved meshes (e.g. each face touches 2 cells or fewer)
- 3D meshes currently only support tetrahedrons

## Application Layer
The application layer is largely unfinished.
Physics models are currently being prototyped to explore and refine the right abstractions.