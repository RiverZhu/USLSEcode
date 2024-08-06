function x = OMP(y, A, P, lambda, on_grid)
% OMP algorithm for argmin_{x} ||Ax - y||^2_2 subject to a P sparse vector x;
% INPUTS:
%   y  - measurements with size M x 1;
%   A  - measurement matrix with size M x N;
%   P  - number of nonzero elements of x;
%   lambda  - dynamic range of ADC;
%   on_grid  - "1" means on_grid; "0" means off_grid (for now, it seems
%   that on-grid OMP performs better);
% OUTPUTs:
%   x  - estimated P sparse vector; 

RD   = @(x,L) round(x./(2*L)).*2*L;

y_r = y;
[~,N] = size(A);
locList = nan(P,1);
gainList = nan(P,1);
normAomp = (sum(abs(A).^2))'; 
for k = 1:P
    Ahyr_energy = (abs(A'*y_r).^2)./normAomp;
    [~,idx] = max(Ahyr_energy);
    locList(k) = idx;
    if on_grid == "0"
        gainList(k) = (A(:,idx))'*y_r/normAomp(idx);
    elseif on_grid == "1"
        gainList(k) = RD((A(:,idx))'*y_r/normAomp(idx), lambda); 
    end
    y_r = y_r - gainList(k)*A(:,idx);
end
x = zeros(N,1);
x(locList) = gainList;
end