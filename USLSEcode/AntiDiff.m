function epsilon = AntiDiff(delta_epsilon, epsilon_true, lambda, V)
% Obtain the anti-difference of delta_epsilon
% INPUTS:
%   delta_epsilon  - the difference of epsilon;
%   epsilon_true  - first element of original epsilon
%   lambda  - dynamic range of ADC
% OUTPUTs:
%   epsilon  - estimated simple function epsilon 

RD = @(x,L) round(x./(2*L)).*2*L; % Round Function

offset_delta = RD(mean(delta_epsilon), lambda);
delta_epsilon_new = delta_epsilon - offset_delta;
if max(max(abs(real(delta_epsilon_new))), max(abs(imag(delta_epsilon_new)))) < V + 1e-2
    delta_epsilon = delta_epsilon_new;
end

epsilon = [0; cumsum(delta_epsilon)]; 
offset = RD(mean(epsilon_true - epsilon), lambda); 
epsilon = epsilon + offset;

end