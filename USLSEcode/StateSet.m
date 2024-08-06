function state_set = StateSet(max_value, lambda)
%% Obtain set of all possible states; 
% INPUTS:
%   max_value  - maximum absolute value of both real and imaginary parts of
%   simple function (divide by 2*lambda);
%   lambda  - dynamic range of ADC;
% OUTPUTs:
%   state_set  - the set of all possible states;

state_set_real = -max_value*2*lambda:2*lambda:max_value*2*lambda;
state_set = state_set_real + 1j*(state_set_real');
state_set = state_set(:);
end