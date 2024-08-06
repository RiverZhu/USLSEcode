function x = DPBackward(x_vec_idx, x_opt, x_set)
% backward process of dynamic programming;
% INPUTS:
%   x_vec_idx  - optimal x_i for all possible x_{i+1}:x_{i+K};
%   x_opt  - optimal x_{i+1}:x_{i+K};
%   x_set  - candidates for each parameter;
% OUTPUTs:
%   x  - optimal x_i; 

K = length(x_opt);
candi_num = length(x_set);
[~, idx_opt] = ismember(x_opt, x_set);
idx = (candi_num.^(0 : K - 1)) * (idx_opt - 1) + 1;
x_idx = x_vec_idx(idx);
x = x_set(x_idx);
end