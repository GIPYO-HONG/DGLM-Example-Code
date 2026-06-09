function [M, M_I, Dr, Ds] = get_matrices_2d_triangle()
%GET_MATRICES_2D_TRIANGLE Reference P1 matrices on the triangle and interval.

M = [2 1 1; 1 2 1; 1 1 2] / 6;
M_I = [2 1; 1 2] / 3;

Dr = [-1 1 0; -1 1 0; -1 1 0] / 2;
Ds = [-1 0 1; -1 0 1; -1 0 1] / 2;

end
