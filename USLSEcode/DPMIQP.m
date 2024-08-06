function x = DPMIQP(y, A, x_set, K, x_true)
%% Dynamic programming for argmin_{x} ||Ax - y||^2_2 subject to x_i in x_set and K-th order Markov assumption;
% INPUTS:
%   y  - measurements with size M x 1;
%   A  - measurement matrix with size M x N;
%   x_set  - candidates for each parameter;
%   K  - K-th order Markov assumption;
%   lambda  - dynamic range of ADC; 
%   x_true - the real value of parameters;
% OUTPUTs:
%   x  - estimated parameters; 

%% common functions;
Loss = @(x,Q,b,mu) real(x' * Q * x - 2 * real(b' * x));

%% initialization;
Q = A' * A; b = A' * y;
[~,N] = size(A);  % number of parameters;
candi_num = length(x_set);  % number of candidates;
DP = zeros(candi_num^K, N - K);  % initialize the DP matrix;
x_idx = zeros(candi_num^K, N - K);  % initialize the optimal index;
x = zeros(N, 1);  % initialize the estimated parameters;
eps = 1e-4;  % a small number;

%% dynamic programming algorithm;
% forward process;
for i = 1 : N - K - 1
    Q_tilde = [Q(i, i : i + K); [Q(i + 1 : i + K, i), zeros(K)]]; 
    [DP(:, i + 1), x_idx(:, i)] = DPForward(Q_tilde, b(i), DP(:, i), x_set);
end
% special handling for the last step;
Q_tilde = Q(N - K : N, N - K : N);
[xend, loss] = DPForwardLast(Q_tilde, b(N - K: N), DP(:, N - K), x_set); 
% backward process;
x(N - K: end) = xend;
for i = N - K - 1 : -1 : 1
    x(i) = DPBackward(x_idx(:, i), x(i + 1: i + K), x_set);
end

%% check if DP is correct;
Q_approx = Q - triu(Q, K + 1) - tril(Q, -(K + 1));
if abs(Loss(x, Q_approx, b) - loss) > eps
    error('Loss of DP algorithm is not consistant!')
end
x_true_range = max(max(abs(real(x_true))), max(abs(imag(x_true))));
x_set_range = max(max(abs(real(x_set))), max(abs(imag(x_set))));
if Loss(x_true, Q_approx, b) - Loss(x, Q_approx, b) < -eps && x_true_range <= x_set_range + eps
    error('DP algorithm does not obtain the optimal solution!')
end

end