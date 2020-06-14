%% SetVariable()
%  Set the variables.
%
%  Syntax
%
%  Descriptions
%
%%
function SetVariable()
load('para.mat')
global SOLNEARZERO

% Set space
rDense = 0:drDense:rend;

% Set the ode solver
odeinput4a.opts = odeset('RelTol', 1e-5, 'MaxStep', MaxStep4a, 'Events', @StopEvents);
odeinput4a.range = [0, rend];
odeinput4a.ic = ica;
odeinput4alpha.opts = odeset('RelTol', 1e-5, 'MaxStep', MaxStep4alpha);
odeinput4alpha.range = [0, rend];

% Extrapolation of phi
A(1, :) = [drDense^2 1];
A(2, :) = [4*drDense^2 1];
SOLNEARZERO = inv(A);

save('vari.mat')