# FEMto
Finite element method simulation written completely from scratch in Odin.

> [!NOTE]
> FEMto is in early development.


## Overview
FEMto is a finite element solver designed to solve partial differential equations representing common physics. Currently, it focuses on
heat conduction, including nonlinear effects such as temperature-dependent materials and radiative boundaries,
but the core library is ready for vector-valued and coupled PDEs.

FEMto is designed to be small and simple, written in a plain procedural style without relying on third-party libraries.
A text-based configuration system wraps a subset of its functionality to allow rapid experimentation and iteration.

> [!NOTE]
> The config interface is not intended to expose every feature of FEMto. To fully leverage its capabilities,
> users can implement custom materials, sources, or boundary conditions directly in code.


## Features
- 1D, 2D, 3D linear or curved geometries from [Gmsh](https://gmsh.info/).
- Visualization via [Paraview](https://www.paraview.org/).
- Steady-state and transient simulations.
- Text based config using a human-friendly JSON superset.
- Support for time, space, and solution dependent material models, source terms and boundary conditions.
- Non linear support via Newton *Picard-style supported through lagged values in the weak form.*
- 1st and 2nd order solution accuracy.

### Mesh Restrictions
- Only tetrahedrons for 3D so far.
- Meshes should be well-behaved.
- Every entity including domain entities needs a physical name.
- **Version 2.2 Binary only**.

### Current Limitations
- Linear solver is limited. We currently rely on conjugate gradient with no preconditioning. Complex problems will struggle.
- Dirichlet boundaries currently do not time/space varying values yet.
- No way to read initial conditions in from a file from the config.
- Limited built-in sources, BCs and materials other than plain constant ones.


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
name = "conduction_sim" //name is just for output files

mesh = "../meshes/cantilever.msh" // relative to config file

model = Conduction // the rest of our config is now specific to a conduction simulation

// sections come from names on mesh entities, by breaking up your domain into multiple sections
// in Gmsh you can have piece-wise constant mat props, sources and initial conditions
sections = {
    beam = { 
        material = {
            type = "constant",
            conductivity = 50000, //unphysical just to treat the body as isothermal
            density = 7800,
            specific_heat = 490,
        }
    }
}

fields = {
    temperature = {
        order = Linear, // solution order
        boundaries = {
            // the name comes from the mesh, could be anything "chicken" = {type = ...}
            "constraint" = {type = "radiative", emissivity = 0.9, ambient_temp = 517},
            "free" = {type = "radiative",emissivity = 0.9, ambient_temp = 517 }
        }
        initial_conditions = {
            // heres the section name again.
            beam = {
                type = "constant",
                temperature = 400,
            }
        }
    }
}


output = {
    directory = "../output",
    frequency = 1, // only meaningful for transient simulations, how often to output results
}

// optional
transient = {
    end = 100 // 100 units
    timestep = 1, // 1 unit at a time
    //start = xxx can also define a start time
}

// optional
linear_solver = {
    rtol = 1e-6
    max_iterations = 5000
}

// optional
non_linear_solver = {
    rtol = 1e-5
}

```

## Contributions
Welcome, especially if you're a wizard with multi-grid solvers.

## License
See LICENSE file.
