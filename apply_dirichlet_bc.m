function [K_global, F_global] = apply_dirichlet_bc( ...
    K_global, F_global, ...
    c4n, ind4e, inddb, g_func)
%APPLY_DIRICHLET_BC Strongly impose Lambda = g on boundary edge DoFs.

NumElements = size(ind4e, 1);
m = 3;
n = 2;
U_DoF_Size = NumElements * m;

vertex_loc = [2 3 1];
local_edges = {[2 3], [3 1], [1 2]};

edges_local = [
    ind4e(:,vertex_loc(1)), ind4e(:,vertex_loc(2));
    ind4e(:,vertex_loc(2)), ind4e(:,vertex_loc(3));
    ind4e(:,vertex_loc(3)), ind4e(:,vertex_loc(1))];

[~, ~, edge_indices] = unique(sort(edges_local, 2), 'rows');
element_to_edge = reshape(edge_indices, NumElements, 3);

is_boundary_node = false(size(c4n, 1), 1);
is_boundary_node(inddb) = true;

num_boundary_edges = 0;

for elem = 1:NumElements
    elem_nodes = ind4e(elem, :);

    for edge = 1:3
        loc = local_edges{edge};
        edge_nodes = elem_nodes(loc);

        if is_boundary_node(edge_nodes(1)) && is_boundary_node(edge_nodes(2))
            num_boundary_edges = num_boundary_edges + 1;
        end
    end
end

boundary_lambda_dofs = zeros(1, n*num_boundary_edges);
boundary_g_vals = zeros(1, n*num_boundary_edges);
count = 0;

for elem = 1:NumElements
    elem_nodes = ind4e(elem, :);

    for edge = 1:3
        loc = local_edges{edge};
        edge_nodes = elem_nodes(loc);

        if ~(is_boundary_node(edge_nodes(1)) && is_boundary_node(edge_nodes(2)))
            continue;
        end

        global_edge = element_to_edge(elem, edge);
        lambda_idx = U_DoF_Size + ((global_edge-1)*n + 1:global_edge*n);

        if edge_nodes(1) > edge_nodes(2)
            lambda_idx = fliplr(lambda_idx);
        end

        vals = [
            g_func(c4n(edge_nodes(1), 1), c4n(edge_nodes(1), 2)), ...
            g_func(c4n(edge_nodes(2), 1), c4n(edge_nodes(2), 2))];

        boundary_lambda_dofs(count+1:count+n) = lambda_idx;
        boundary_g_vals(count+1:count+n) = vals;
        count = count + n;
    end
end

[boundary_lambda_dofs, ia] = unique(boundary_lambda_dofs);
boundary_g_vals = boundary_g_vals(ia);

g_vec = sparse(boundary_lambda_dofs, ones(size(boundary_lambda_dofs)), ...
               boundary_g_vals', size(K_global, 2), 1);
F_global = F_global - K_global * g_vec;

for idx = 1:length(boundary_lambda_dofs)
    row = boundary_lambda_dofs(idx);
    K_global(row, :) = 0;
    K_global(:, row) = 0;
    K_global(row, row) = 1;
    F_global(row) = boundary_g_vals(idx);
end

end
