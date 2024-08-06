function G = Cartesian(args)
%% Obtain the cartesian product of several sets;
% INPUTS:
%   args  - a cell contains several sets
% OUTPUTS:
%   G  - a matrix to represent the cartesian product

n = length(args);
[F{1:n}] = ndgrid(args{:});
for i=n:-1:1
    G(:,i) = F{i}(:);
end
end