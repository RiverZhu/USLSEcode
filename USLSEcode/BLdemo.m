% Demo for the USLSE algorithm for bandlimited signals via unlimited sampling
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

K = @(x) sinc(x);  % Bandlimited Kernel

%% parameters setting
SNR = 25;   % signal_to_noise ratio
gamma = 10; % oversampling factor
lambda = 0.5;   % dynamic range of ADC

%% generate the signal
t  = -20:0.005:20 - 0.005; t = t(:);
Z = 10 * (rand(numel(t), 1) + 1i .* rand(numel(t), 1));

F = fft(K(t)) .* Z;
sig_I = real(ifft(F));
F_any = hilbert(sig_I); 
F_any_FFT = fft(F_any);
ts = 1 / gamma;
sam_int = round(ts / 0.005);
locs = 1:sam_int:numel(t);

x = F_any(locs);    % generate the noiseless signal 
N = length(x);
sigma = norm(x)^2 / N / (10^(SNR / 10));
noise = sqrt(sigma/2) .* (randn(N, 1) + 1j * randn(N, 1));  % generate the noise
g = x + noise; % unfolded samples
y = Mod(g, lambda); % modulo measurements
epsilon = g - y;  % simple function
disp(['The maximum folding times is ', num2str(ceil((CMAX(g) - lambda) / 2 / lambda)), '.']);

%% Hyperparameters of algorithms
DP_hyper = struct('beta', 0.02, 'P', 4, 'V', 1);  % hyperparameters of DP algorithm

%% DP algorithm
F = dftmtx(N-1) / sqrt(N-1);  % DFT matrix
del_epsilon = CD(epsilon); del_y = CD(y); del_g = CD(g); del_x = CD(x);
fre_del_epsilon = F * del_epsilon; fre_del_y = F * del_y; fre_del_g = F * del_g; fre_del_x = F * del_x;
M_set = floor((N-1) * (1 / gamma + DP_hyper.beta)) + 2 : floor((N - 1) * (1 - DP_hyper.beta)) + 1;
del_epsilon_set = StateSet(DP_hyper.V, lambda);

del_epsilon_est = zeros(N-1, 1);
for iter = 1 : 3
    del_epsilon_add_est = DPMIQP((-fre_del_y(M_set) - F(M_set, :) * del_epsilon_est), F(M_set, :), del_epsilon_set, DP_hyper.P, del_epsilon);
    del_epsilon_est = del_epsilon_est + del_epsilon_add_est;
    del_epsilon_add_est = OMP((-fre_del_y(M_set) - F(M_set, :) * del_epsilon_est), F(M_set, :), N - 1, lambda, '1');
    del_epsilon_est = del_epsilon_est + del_epsilon_add_est;
end
epsilon = AntiDiff(del_epsilon_est, epsilon, lambda, DP_hyper.V);
DP_g = y + epsilon;

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





