function plot_dglm_solution(sol, c4n, ind4e)
%PLOT_DGLM_SOLUTION Plot the discontinuous P1 solution.

NumElements = size(ind4e, 1);
m = 3;
U_sol = sol(1:NumElements*m);

vertex_loc = [2 3 1];

figure('Color', 'w');
hold on;

for elem = 1:NumElements
    U_idx = ((elem-1)*m + 1):(elem*m);
    elem_nodes = ind4e(elem, :);
    vertices = elem_nodes(vertex_loc);
    vertex_u = U_sol(U_idx(vertex_loc));

    patch( ...
        'XData', c4n(vertices, 1), ...
        'YData', c4n(vertices, 2), ...
        'ZData', vertex_u, ...
        'CData', vertex_u, ...
        'FaceColor', 'interp', ...
        'EdgeColor', [0.3 0.3 0.3], ...
        'LineWidth', 0.5);
end

view(3);
grid on;
box on;
colormap(jet);
colorbar;

xlabel('x', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 11, 'FontWeight', 'bold');
zlabel('u_h(x,y)', 'FontSize', 11, 'FontWeight', 'bold');
title('DGLM solution', 'FontSize', 13, 'FontWeight', 'bold');

rotate3d on;
hold off;

end
