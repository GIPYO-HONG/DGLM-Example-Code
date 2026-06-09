# DGLM Example Code

This repository contains a simple MATLAB implementation of the **Discontinuous Galerkin method with Lagrange Multipliers (DGLM)** for the Poisson equation with non-homogeneous Dirichlet boundary conditions.

The purpose of this repository is to provide a **minimal educational example** of DGLM rather than a complete solver. Only a simple P1 triangular discretization on a uniform mesh is considered.

---

## Problem Description

We solve

$$
-\Delta u = f
\qquad \text{in } \Omega,
$$

with Dirichlet boundary conditions

$$
u = g
\qquad \text{on } \partial \Omega,
$$

on the unit square

$$
\Omega=(0,1)^2.
$$

In the provided example,

$$
f(x,y)=2\pi^2\sin(\pi x)\sin(\pi y),
$$

and

$$
g(x,y)=
\begin{cases}
\sin(x), & y=0,\\
\sin(y), & x=0,\\
0, & x=1 \text{ or } y=1.
\end{cases}
$$

The Dirichlet boundary condition is imposed through edge-based Lagrange multiplier unknowns.

---

## Stabilization Parameter

The stabilization parameter `tau` controls the coupling between neighboring elements.

In this example,

```matlab
h = 1/N;

tau_int = 1/h^2;
tau_bnd = 1/h^2;
```

where

- `tau_int` is used on interior edges,
- `tau_bnd` is used on boundary edges.

The choice

```matlab
tau = 1/h^2
```

was found experimentally to provide a good balance between continuity and conditioning of the linear system. For the test problem considered here, this choice produced relatively small jumps across element interfaces while avoiding unnecessarily large condition numbers.

Although larger values of `tau` generally enforce stronger coupling between neighboring elements, they may also increase the condition number of the system matrix. For this reason, `tau = 1/h^2` was adopted as a reasonable practical choice in this example.

For different source terms, boundary conditions, or meshes, users are encouraged to experiment with different values of `tau` and determine suitable choices for their particular problems.

---

## Usage

All parameters can be modified directly in `main.m`.

### Mesh Resolution

```matlab
N = 8;

Mx = N;
My = N;
```

Increasing `N` generates a finer mesh.

### Stabilization Parameter

```matlab
h = 1/N;

tau_int = 1/h^2;
tau_bnd = 1/h^2;
```

### Run

```matlab
main
```

The code will

1. Generate a triangular mesh.
2. Assemble the DGLM system.
3. Apply Dirichlet boundary conditions.
4. Solve the linear system.
5. Visualize the discontinuous solution.

---

## Repository Structure

```text
main.m                      Main driver script
mesh_fem_2d_triangle.m      Mesh generation
assemble_global_system.m    DGLM system assembly
apply_dirichlet_bc.m        Dirichlet boundary treatment
get_matrices_2d_triangle.m  Reference element matrices
plot_dglm_solution.m        Solution visualization
```

---

## Notes

This repository is intended as a compact example demonstrating

- Poisson equation discretization,
- discontinuous P1 finite elements,
- edge-based Lagrange multipliers,
- stabilization through `tau`,
- and basic solution visualization.

Users are encouraged to modify the source term `f`, boundary condition `g`, mesh resolution, and stabilization parameters to explore the behavior of the DGLM method.
