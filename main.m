clear; clc;

%% Parameters
k = 1;
N = 8;
h = 1 / N;
Mx = N;
My = N;
xl = 0; xr = 1;
yl = 0; yr = 1;

tau_int = 2;
tau_bnd = 2;

%% PDE: -Delta u = f, u = g on boundary
f_func = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);

g_func_ = @(x,y) (y==0).*sin(x) + ... 
           (y==1).*0 + ...           
           (x==0).*sin(y) + ... 
           (x==1).*0;  

%% Mesh
[c4n, ind4e, inddb] = ...
    mesh_fem_2d_triangle(xl, xr, yl, yr, Mx, My, k);

%% Assembly
[K_global, F_global] = assemble_global_system( ...
    c4n, ind4e, ...
    tau_int, tau_bnd, f_func);

%% Dirichlet boundary condition: Lambda = g on boundary edges
[K_global, F_global] = apply_dirichlet_bc( ...
    K_global, F_global, ...
    c4n, ind4e, inddb, g_func_);

%% Solve
sol = K_global \ F_global;

%% Visualization
plot_dglm_solution(sol, c4n, ind4e);
