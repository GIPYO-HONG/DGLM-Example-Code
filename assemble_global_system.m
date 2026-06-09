function [K_global, F_global] = assemble_global_system( ...
    c4n, ind4e, tau_int, tau_bnd, f_func)
%ASSEMBLE_GLOBAL_SYSTEM Assemble the P1 DGLM system for Poisson equation.
%
% Unknown ordering:
%   [elementwise U values; edgewise Lambda values]
%
% The matrix has the block form from the lecture notes:
%   [ A  -B ] [U     ] = [R]
%   [ C  2D ] [Lambda]   [0]

if size(ind4e, 2) ~= 3
    error('assemble_global_system supports P1 elements only.');
end

[M, M_I, Dr, Ds] = get_matrices_2d_triangle();

NumElements = size(ind4e, 1);
m = 3;
n = 2;

vertex_loc = [2 3 1];
local_edges = {[2 3], [3 1], [1 2]};

edges_local = [
    ind4e(:,vertex_loc(1)), ind4e(:,vertex_loc(2));
    ind4e(:,vertex_loc(2)), ind4e(:,vertex_loc(3));
    ind4e(:,vertex_loc(3)), ind4e(:,vertex_loc(1))];

[~, ~, edge_indices] = unique(sort(edges_local, 2), 'rows');
element_to_edge = reshape(edge_indices, NumElements, 3);
NumEdges = max(edge_indices);

edge2elem = zeros(NumEdges, 2);
for elem = 1:NumElements
    for edge = 1:3
        global_edge = element_to_edge(elem, edge);
        if edge2elem(global_edge, 1) == 0
            edge2elem(global_edge, 1) = elem;
        else
            edge2elem(global_edge, 2) = elem;
        end
    end
end
is_boundary_edge = edge2elem(:, 2) == 0;

Total_U_DoF = NumElements * m;
Total_Lambda_DoF = NumEdges * n;

I_A = zeros(NumElements*m*m, 1);
J_A = zeros(NumElements*m*m, 1);
V_A = zeros(NumElements*m*m, 1);

I_B = zeros(NumElements*3*m*n, 1);
J_B = zeros(NumElements*3*m*n, 1);
V_B = zeros(NumElements*3*m*n, 1);

I_C = zeros(NumElements*3*m*n, 1);
J_C = zeros(NumElements*3*m*n, 1);
V_C = zeros(NumElements*3*m*n, 1);

I_D = zeros(NumEdges*n*n, 1);
J_D = zeros(NumEdges*n*n, 1);
V_D = zeros(NumEdges*n*n, 1);

count_A = 0;
count_BC = 0;
count_D = 0;
R_global = zeros(Total_U_DoF, 1);

for elem = 1:NumElements
    U_idx = ((elem-1)*m + 1):(elem*m);
    elem_nodes = ind4e(elem, :);

    xv = c4n(elem_nodes(vertex_loc), 1);
    yv = c4n(elem_nodes(vertex_loc), 2);
    x_all = c4n(elem_nodes, 1);
    y_all = c4n(elem_nodes, 2);

    xr = (xv(1) - xv(3)) / 2;
    xs = (xv(2) - xv(3)) / 2;
    yr = (yv(1) - yv(3)) / 2;
    ys = (yv(2) - yv(3)) / 2;
    detJ = xr * ys - xs * yr;
    absJ = abs(detJ);

    rx =  ys / detJ;
    sx = -yr / detJ;
    ry = -xs / detJ;
    sy =  xr / detJ;

    Dx = rx * Dr + sx * Ds;
    Dy = ry * Dr + sy * Ds;
    A_local = absJ * (Dx' * M * Dx + Dy' * M * Dy);

    B_edges = cell(3, 1);
    C_edges = cell(3, 1);
    D_edges = cell(3, 1);

    for edge = 1:3
        global_edge = element_to_edge(elem, edge);
        tau = tau_int;
        if is_boundary_edge(global_edge)
            tau = tau_bnd;
        end

        loc = local_edges{edge};
        edge_nodes = elem_nodes(loc);
        p1 = c4n(edge_nodes(1), :);
        p2 = c4n(edge_nodes(2), :);
        edge_vec = p2 - p1;
        edge_len = norm(edge_vec);
        normal = [edge_vec(2), -edge_vec(1)] / edge_len;

        normal_derivative = Dx(loc, :) * normal(1) + Dy(loc, :) * normal(2);
        P = (edge_len / 2) * normal_derivative' * M_I;
        D = (tau / 2) * M_I;

        P_bar = zeros(m, m);
        D_bar = zeros(m, m);
        P_bar(:, loc) = P;
        D_bar(loc, loc) = D;

        A_local = A_local - P_bar' + P_bar + D_bar;
        B_edges{edge} = P + D_bar(:, loc);
        C_edges{edge} = P' - D_bar(loc, :);
        D_edges{edge} = D;
    end

    f_vals = zeros(m, 1);
    for i = 1:m
        f_vals(i) = f_func(x_all(i), y_all(i));
    end
    R_global(U_idx) = R_global(U_idx) + absJ * M * f_vals;

    for i = 1:m
        for j = 1:m
            count_A = count_A + 1;
            I_A(count_A) = U_idx(i);
            J_A(count_A) = U_idx(j);
            V_A(count_A) = A_local(i, j);
        end
    end

    for edge = 1:3
        global_edge = element_to_edge(elem, edge);
        lambda_idx = ((global_edge-1)*n + 1):(global_edge*n);
        loc = local_edges{edge};
        edge_nodes = elem_nodes(loc);

        if edge_nodes(1) > edge_nodes(2)
            lambda_idx = fliplr(lambda_idx);
        end

        B = B_edges{edge};
        C = C_edges{edge};
        D = D_edges{edge};

        for i = 1:m
            for j = 1:n
                count_BC = count_BC + 1;
                I_B(count_BC) = U_idx(i);
                J_B(count_BC) = lambda_idx(j);
                V_B(count_BC) = -B(i, j);

                I_C(count_BC) = lambda_idx(j);
                J_C(count_BC) = U_idx(i);
                V_C(count_BC) = C(j, i);
            end
        end

        if edge2elem(global_edge, 1) == elem
            for i = 1:n
                for j = 1:n
                    count_D = count_D + 1;
                    I_D(count_D) = lambda_idx(i);
                    J_D(count_D) = lambda_idx(j);
                    V_D(count_D) = 2 * D(i, j);
                end
            end
        end
    end
end

A_global = sparse(I_A(1:count_A), J_A(1:count_A), V_A(1:count_A), ...
                  Total_U_DoF, Total_U_DoF);
B_global = sparse(I_B(1:count_BC), J_B(1:count_BC), V_B(1:count_BC), ...
                  Total_U_DoF, Total_Lambda_DoF);
C_global = sparse(I_C(1:count_BC), J_C(1:count_BC), V_C(1:count_BC), ...
                  Total_Lambda_DoF, Total_U_DoF);
D_global = sparse(I_D(1:count_D), J_D(1:count_D), V_D(1:count_D), ...
                  Total_Lambda_DoF, Total_Lambda_DoF);

K_global = [A_global, B_global; C_global, D_global];
F_global = [R_global; zeros(Total_Lambda_DoF, 1)];

end
