function [x, loss] = DPForwardLast(Q_tilde, b, DP, x_set)
%% forward process of dynamic programming of the last step;
% INPUTS:
%   Q_tilde  - coefficients of quadratic terms;
%   b  - coefficients of linear terms;
%   DP  - solution of DP for the last time step;
%   x_set  - candidates for each parameter;
%   lambda  - dynamic range of ADC;
% OUTPUTs:
%   x  - estimates of last K parameters;
%   loss  - final loss;

[~, R] = size(Q_tilde);
candi_num = length(x_set);
x_cell = cell(1, R); x_cell(:) = {x_set};
x_vec = Cartesian(x_cell);
DP_rep = repmat(DP, candi_num, 1);
loss_vec = real(sum((conj(x_vec) * Q_tilde) .* x_vec, 2)) - 2 * real(conj(x_vec) * b) + DP_rep;
[loss, idx] = min(loss_vec);
x = x_vec(idx, :);
end