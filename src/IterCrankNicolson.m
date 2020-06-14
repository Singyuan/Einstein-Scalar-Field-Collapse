%% IterCrankNicolson
%  Solve the next time step scalar field by Iterative Crank Nicolson
%  method.
%  Please refer to Alcubierre et al. (2000) and Teukolsky (2000)
%
%  Syntax
%
%  Descriptions
%
%%
function [U, V, W] = IterCrankNicolson(U, V, W, v, x, dx, dt)
% Declare the global variables
global SOLNEARZERO ICN

% Declare the persistent variables
% The variables Udiffu and Vdiffu are diffusion terms.
% The variables dUold dVold are previous derivative.
persistent Udiffu Vdiffu dUold dVold

% The first sub-step approximate the next time step
if ICN == 1
    % Prepare the diffusion terms
    Udiffu = myfourdudr(U(1, :));
    Vdiffu = myfourdudr(V(1, :));
    
    % Save origin time step differentiation of phi and pi
    % Correct inner boundary by specific method
    dUold = my5orderdudr_Gust(x.^2.*v.*U(1, :), dx);
    dUold = (1./x.^2).*dUold;
    dUold(1:5) = Evansinnerbdy(x, v, U); % Correct inner boundary
    dVold = my5orderdudr_Gust(v.*V(1, :), dx);
    
    % Approxiate phi and pi at the next step.
    % Add the diffusion term.
    V(2, :) = V(1, :)+dt*dUold+1/32*Vdiffu;
    U(2, :) = U(1, :)+dt*dVold+1/32*Udiffu;
    
    % Correct the inner boundary
    % pi is symmetric to origin so I use two points to meet 4 order.
    temp = (SOLNEARZERO*V(2, 2:3)');
    V(2, 1) = temp(2);
    U(2, 1) = 0;
    
    % Approximate scalar field at the next step.
    % Solve sf by Sommefeld condition
    W(2, :) = W(1, :)+dt*v.*V(1, :);
    W(2, end-1:end) = SommerfeldOuterBdy(x(end-1:end), W(:, end-3:end), dt, dx);
    
    % Correct the outer boundary of phi
    dW = my5orderdudr_Gust(W(2, :), dx);
    U(2, end-2:end) = dW(end-2:end);
    
% The second and third sub-step of correct to next step
else
    % Compute the next step differentiation
    % Correct inner boundary by specific method
    dU = my5orderdudr_Gust(x.^2.*v.*U(2, :), dx);
    dU = (1./x.^2).*dU;
    dUold(1:5) = Evansinnerbdy(x, v, U); % Correct inner boundary
    dV = my5orderdudr_Gust(v.*V(2, :), dx);
    
    % "Sneak" phi and pi at the next step.
    V(2, :) = V(1, :)+0.5*dt*(dU+dUold);
    U(2, :) = U(1, :)+0.5*dt*(dV+dVold);
    
    % Correct the inner boundary
    temp = (SOLNEARZERO*V(2, 2:3)');
    V(2, 1) = temp(2);
    U(2, 1) = 0;
    
    % Approxiate sf at the next step.
    % Solve sf by Sommefeld condition
    W(2, :) = W(1, :)+0.5*dt*v.*(V(1, :)+V(2, :));
    W(2, end-1:end) = SommerfeldOuterBdy(x(end-1:end), W(:, end-3:end), dt, dx);
    
    % Correct the outer boundary of phi
    dW = my5orderdudr_Gust(W(2, :), dx);
    U(2, end-2:end) = dW(end-2:end);
end
    
if ICN == 3   
    % Update
    U(1, :) = U(2, :);
    V(1, :) = V(2, :);
    W(1, :) = W(2, :);
end
end