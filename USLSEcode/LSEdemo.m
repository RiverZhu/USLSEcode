% Demo for the USLSE algorithm for line spectral signals via unlimited sampling
% For details, refer to the paper
% Q. Zhang, J. Zhu*, F. Qu and D. W. Soh, Line Spectral Estimation via Unlimited Sampling,
% IEEE Transactions on Aerospace and Electronic Systems, 2024 
% This code is written by Qi Zhang and Jiang Zhu. If you have any problems, please feel free to contact
% jiangzhu16@zju.edu.cn.

clear variables;
rng(123)

%% common functions 
Frac = @(x)  x - floor(x);  % Fractional Part

CD = @(x) diff(x, 1);  % difference operator

RD = @(x,L) round(x ./ (2 * L)) .* 2 * L;  % Round Function

Mod = @(f,T) 2 * T .* (Frac(real(f) ./ (2 .* T) + 0.5) -0.5 ... 
    + 1j * (Frac(imag(f) ./ (2 .* T) + 0.5) - 0.5) );  % Origin Centered Modulo Function

CMAX = @(x) max(max(abs(real(x(:)))), max(abs(imag(x(:)))));  % Obtain the maximum value of the absolute value of the real part and imaginary part of x;


%% parameters setting
N = 512;    % number of samples
K = 3;      % number of signals
omega = [0.31; 0.42; 0.5];    % frequencies
SNR = 25;   % signal_to_noise ratio
gamma = 10; % oversampling factor
lambda = 0.5;   % dynamic range of ADC
sigma_c = 2;    % variance of the complex weights


%% generate the signal
c = sqrt(sigma_c / 2) * randn(K, 1) + sqrt(sigma_c / 2) * 1j * randn(K, 1);  
x = exp(1j * (0: N-1)' * omega') * c;   % generate the noiseless signal 
sigma = norm(x, 'fro')^2 / N / (10^(SNR / 10)); 
noise = sqrt(sigma / 2) .* (randn(N, 1) + 1j * randn(N, 1));  % generate the noise
g = x + noise;  % unfolded samples
y = Mod(g, lambda); % modulo measurements
epsilon = g - y;  % simple function
disp(['The maximum folding times is ', num2str(ceil((CMAX(g) - lambda) / 2 / lambda)), '.']);

%% Hyperparameters of algorithms
DP_hyper = struct('beta', 0.02, 'P', 4, 'V', 1);  % hyperparameters of DP algorithm

%% DP algorithm
F = dftmtx(N-1)/sqrt(N-1);  % DFT matrix
del_epsilon = CD(epsilon); del_y = CD(y); del_g = CD(g); del_x = CD(x);
fre_del_epsilon = F * del_epsilon; fre_del_y = F * del_y; fre_del_g = F * del_g; fre_del_x = F * del_x;
M_set = floor((N-1)*(1/gamma + DP_hyper.beta)) + 2 : floor((N-1)*(1-DP_hyper.beta)) + 1;
del_epsilon_set = StateSet(DP_hyper.V, lambda);

del_epsilon_est = zeros(N-1, 1);
for iter = 1 : 3
    del_epsilon_add_est = DPMIQP((-fre_del_y(M_set) - F(M_set, :) * del_epsilon_est), F(M_set, :), del_epsilon_set, DP_hyper.P, del_epsilon);
    del_epsilon_est = del_epsilon_est + del_epsilon_add_est;
    del_epsilon_add_est = OMP((-fre_del_y(M_set)- F(M_set, :) * del_epsilon_est), F(M_set, :), N-1, lambda, '1');
    del_epsilon_est = del_epsilon_est + del_epsilon_add_est;
end
epsilon = AntiDiff(del_epsilon_est, epsilon, lambda, DP_hyper.V);
DP_g = y + epsilon;

%% Perform LSE
[omega_nomp, gain_nomp, ~] = KMNOMP(DP_g, eye(N), 0.01, K);  % use NOMP to solve LSE

%% Figures
set(groot, 'defaultAxesFontSize', 22);
f1 = figure(1);
set(gcf, 'position', [0 0 1000 800]);

subplot(2,1,1)
box on
hold on
plot(1:N, real(g), '-', 'LineWidth', 2, 'Color', [0.855, 0.647, 0.125]);
stem(1:N, real(y),'-','LineWidth', 1,  'MarkerSize', 4, 'Color', [0.3, 0.3, 0.3]);
plot(1:N, real(DP_g), '-.k', 'LineWidth', 1.5);
plot([1, N], [lambda, lambda], '-.k', 'LineWidth', 1);
plot([1, N], [-lambda, -lambda], '-.k', 'LineWidth', 1);
plot([1, N], [0, 0], '-k', 'LineWidth', 1);
fill([1:N, fliplr(1:N)], [real(g)', transpose(zeros(size(real(g))))], [0.855, 0.647, 0.125], 'FaceAlpha', 0.2, 'EdgeColor','none');     
legend('Ground Truth', 'Modulo Samples', 'USLSE', 'Location', 'NorthEast');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 N])
ylim([min(real(g)) - 0.1, max(real(g)) + 0.1])
title('Real Part')

subplot(2,1,2)
box on
hold on
plot(1:N, imag(g), '-', 'LineWidth', 2, 'Color', [0.855, 0.647, 0.125]);
stem(1:N, imag(y), '-', 'LineWidth', 1,  'MarkerSize', 4, 'Color', [0.3, 0.3, 0.3]);
plot(1:N, imag(DP_g), '-.k', 'LineWidth', 1.5);
plot([1, N], [lambda, lambda], '-.k', 'LineWidth', 1);
plot([1, N], [-lambda, -lambda], '-.k', 'LineWidth', 1);
plot([1, N], [0, 0], '-k', 'LineWidth', 1);
fill([1:N, fliplr(1:N)], [imag(g)', transpose(zeros(size(imag(g))))], [0.855, 0.647, 0.125], 'FaceAlpha', 0.2, 'EdgeColor','none');  
xlabel('Samples');
ylabel('Amplitude');
xlim([0 N])
ylim([min(imag(g)) - 0.1, max(imag(g)) + 0.1])
title('Imaginary Part')

f2 = figure(2);
set(gcf, 'position', [0 0 500 400]);
polarplot(omega, abs(c), 'bo', 'MarkerSize', 12, 'LineWidth', 1);
hold on
polarplot(omega_nomp, abs(gain_nomp)/sqrt(N), 'r*', 'MarkerSize', 12, 'LineWidth', 1);
legend('Truth','USLSE','Fontsize',18);
set(gca, 'FontSize', 18);








