# FEMto
Finite element method simulation written completely from scratch in Odin.

> [!NOTE]
> FEMto is in early development.


## Overview
FEMto is a finite element solver designed to solve partial differential equations representing common physics. Currently, it focuses on
heat conduction and small-strain elasticity. With support for non-linear, time, and space dependant properties.

FEMto is designed to be small and simple, written in a plain procedural style without relying on third-party libraries.
A text-based configuration system wraps a subset of its functionality to allow rapid experimentation and iteration.

> [!NOTE]
> The config interface is not intended to expose every feature of FEMto. To fully leverage its capabilities,
> users can implement custom materials, sources, boundary conditions or even new PDEs directly in code.


## Features
- 1D, 2D, 3D linear or curved geometries from [Gmsh](https://gmsh.info/).
- Visualization via [Paraview](https://www.paraview.org/).
- Steady-state and transient simulations.
- Text based config using a human-friendly JSON superset.
- Support for time, space, and solution dependent material models, source terms and boundary conditions.
- Non linear support via Newton *Picard-style supported through lagged values in the weak form.*
- 1st and 2nd order solution accuracy.

### Mesh Restrictions
- Meshes should be well-behaved.
- Every entity including domain entities needs a physical name.
- **Version 2.2 Binary only**.

### Current Limitations
- Linear solver is limited. We currently rely on bi-conjugate gradient stabilized with no preconditioner.
- Dirichlet boundaries currently do not support time/space varying values yet.
- No way to read initial conditions in from a file from the config.
- Limited built-in sources and materials other than plain constant and linear ones.

## Validation
Nothing yet, mostly qualitative checks. Elasticity matches euler-bernoulli theory for cantilever beams.


## Quickstart

```sh
git clone https://github.com/Rwn-A/FEMto femto

cd femto

odin build src/app -o:speed -out:femto.exe

./femto.exe ./configs/conduction.mjson
```

## Architecture
The project is split into 3 distinct parts.
- **Core:** The base FEM and linear algebra code that supports the physics.
Other than `src/fe_core/layout.odin`, all of this code is method-agnostic CG, DG, HDG, etc. are all viable given the appropriate implementation.
- **Models:** The actual physics layer. Exposes an interface for implementing new weak forms, materials, boundaries and sources.
- **App:** Responsible for reading config and constructing a model. Not intended to expose every feature, for full control, implement directly against the Models interface.

## Configuration
Configuration is unstable while development is ongoing. At the time of writing the below example should be a good place to start.
```
name = "part_elastic"
mesh = "../meshes/cantilever_structured.msh"
model = Elasticity
sections = {
    beam = {
        material = {
            type = "constant",
            elastic_modulus = 2e11,
            poissons_ratio = 0.3,
        }
        sources = [
        ]
    }
}
fields = {
    displacement = {
        order = Linear,
        boundaries = {
            "end_left" = [{type = "fixed", displacement=[0, 0, 0]}],
            "end_right" = [{type = "traction", traction = [0, 0, 2e6]}]
        }
    }
}

linear_solver = {
    rtol = 1e-4
    max_iterations = 5000
}

non_linear_solver = {
    rtol = 1e-3
}

output = {
    directory = "../output",
    frequency = 1,
}
```

## Contributions
Welcome, especially if you're a wizard with multi-grid solvers.

## License
See LICENSE file.
