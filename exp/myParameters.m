function myParameters()
% set the parameter of constant

global SOLNEARZERO
fg = 0.8; % CFL constant
T = 10; % total time


rend = 30; % end point of r
drDense = 0.1; % gird size
rDense = 0:drDense:rend;



% initial scalar field for Marse
phi0 = 0.001;
pow = 2;
r0 = 5;
delta = 1;
q = 2;

% C.phi0 = 0.01;
% C.pow = 0;
% C.r0 = 20;
% C.delta = 5;
% C.q = 2;


% solve ode
MaxStep4a = 0.2;
MaxStep4alpha = 0.2;
ica = 1; % initial value of a

A(1, :) = [drDense^2 1];
A(2, :) = [4*drDense^2 1];
SOLNEARZERO = inv(A);

save('para')
end