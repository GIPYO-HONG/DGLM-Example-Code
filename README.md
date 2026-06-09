# DGLM Example Code

This repository provides a simple MATLAB implementation of the **Discontinuous Galerkin method with Lagrange Multipliers (DGLM)** for the Poisson equation with non-homogeneous Dirichlet boundary conditions.

## Problem Description

We solve the Poisson equation

\[
-\Delta u = f \qquad \text{in } \Omega,
\]

with Dirichlet boundary conditions

\[
u = g \qquad \text{on } \partial\Omega,
\]

on the unit square domain

\[
\Omega=(0,1)^2.
\]

In the current example, the source term is chosen as

\[
f(x,y)=2\pi^2\sin(\pi x)\sin(\pi y),
\]

and the boundary condition is prescribed by

\[
g(x,y)=
\begin{cases}
\sin(x), & y=0,\\
\sin(y), & x=0,\\
0, & x=1 \text{ or } y=1.
\end{cases}
\]

The Dirichlet boundary condition is imposed through a Lagrange multiplier defined on mesh edges.

---

## DGLM Formulation

The solution is approximated using **discontinuous P1 finite elements**, meaning that the numerical solution is allowed to be discontinuous across element interfaces.

To weakly enforce continuity and boundary conditions, an additional unknown (the Lagrange multiplier \(\lambda\)) is introduced on element edges.

The resulting global system has the block structure

\[
\begin{bmatrix}
A & -B \\
C & 2D
\end{bmatrix}
\begin{bmatrix}
U \\
\Lambda
\end{bmatrix}
=
\begin{bmatrix}
R \\
0
\end{bmatrix},
\]

where

- \(U\) denotes the element-wise discontinuous solution,
- \(\Lambda\) denotes the edge-based Lagrange multiplier.

---

## Stabilization Parameter

The parameter \(\tau\) controls the stabilization of the DGLM formulation.

- Larger values of \(\tau\) impose stronger coupling across interfaces.
- Smaller values lead to weaker stabilization and may affect robustness.
- Proper choices of \(\tau\) can improve stability and accuracy.

The implementation allows different stabilization parameters for

- interior edges (`tau_int`)
- boundary edges (`tau_bnd`)

---

## Usage

All simulation parameters can be modified directly in `main.m`.

### Mesh Resolution

```matlab
N = 8;

Mx = N;
My = N;
```

This generates a uniform triangular mesh on the unit square.

Increasing `N` produces a finer mesh and generally improves the approximation quality.

### Stabilization Parameters

```matlab
tau_int = 2;
tau_bnd = 2;
```

- `tau_int` : stabilization parameter for interior edges.
- `tau_bnd` : stabilization parameter for boundary edges.

### Run

Execute

```matlab
main
```

The code will

1. Generate the mesh.
2. Assemble the DGLM system.
3. Impose Dirichlet boundary conditions.
4. Solve the resulting linear system.
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

This code is intended as a minimal educational example demonstrating

- Discontinuous Galerkin discretization,
- Lagrange multiplier enforcement of boundary conditions,
- Edge-based hybrid unknowns,
- Stabilization through the parameter \(\tau\),

for the Poisson equation on triangular meshes.
