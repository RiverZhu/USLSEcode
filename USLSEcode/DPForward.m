function [DP_next, x_idx] = DPForward(Q_tilde, b, DP, x_set)
% forward process of dynamic programming;
% INPUTS:
%   Q_tilde  - coefficients of quadratic terms;
%   b  - coefficients of linear terms;
%   DP  - solution of DP for the last time step;
%   x_set  - candidates for each parameter;
% OUTPUTs:
%   DP_next  - solution of DP for the next time step;
%   x_idx  - the record of optimal state;

[~, R] = size(Q_tilde);
candi_num = length(x_set);
x_cell = cell(1, R); x_cell(:) = {x_set};
x_vec = Cartesian(x_cell);
DP_rep = repmat(DP, candi_num, 1);
loss = real(sum((conj(x_vec) * Q_tilde) .* x_vec, 2)) - 2 * real(conj(b) * x_vec(:, 1)) + DP_rep;
loss = reshape(loss, candi_num, []); loss = loss'; 
[DP_next, x_idx] = min(loss, [], 2);
end