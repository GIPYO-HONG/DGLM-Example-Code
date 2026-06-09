function [AE, BE, CE, DE] = build_matrices(M, M_I, Dr, Ds, tau)
%BUILD_MATRICES Reference P1 DGLM matrices from the lecture-note formulas.

m = 3;
edge_nodes = {[2 3], [3 1], [1 2]};
ref_nodes = [
    -1 -1;
     1 -1;
    -1  1];

S = Dr' * M * Dr + Ds' * M * Ds;

AE = S;
BE = zeros(m, m);
CE = zeros(m, m);
DE = zeros(m, m);

for edge = 1:3
    loc = edge_nodes{edge};
    p1 = ref_nodes(loc(1), :);
    p2 = ref_nodes(loc(2), :);
    edge_vec = p2 - p1;
    edge_len = norm(edge_vec);
    normal = [edge_vec(2), -edge_vec(1)] / edge_len;

    normal_derivative = Dr(loc, :) * normal(1) + Ds(loc, :) * normal(2);
    P = (edge_len / 2) * normal_derivative' * M_I;
    D = (tau(edge) / 2) * M_I;

    P_bar = zeros(m, m);
    D_bar = zeros(m, m);
    P_bar(:, loc) = P;
    D_bar(loc, loc) = D;

    AE = AE - P_bar' + P_bar + D_bar;
    BE = BE + P_bar + D_bar;
    CE = CE + P_bar' - D_bar;
    DE = DE + D_bar;
end

end
