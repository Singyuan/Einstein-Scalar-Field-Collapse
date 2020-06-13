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
% The variables dUold dVold and dbdyWold are previous derivative.
% The variable dbdyWold is derivative boundary of scalar field.
persistent Udiffu Vdiffu dUold dVold dbdyWold

% The first sub-step approximate the next time step
if ICN == 1
    % Prepare the diffusion terms
    Udiffu = myfourdudr(U(1, :));
    Vdiffu = myfourdudr(V(1, :));
    
    % Save origin time step differentiation of phi and pi
    dUold = my5orderdudr_Gust(x.^2.*v.*U(1, :), dx);
    dUold = (1./x.^2).*dUold;
    dUold(1) = 3*(x(2)^2*v(2)*U(1, 2)-x(1)^2*v(1)*U(1, 1))/(x(2)^3-x(1)^3);
    dUold(2:5) = 3*(x(3:6).^2.*v(3:6).*U(1, 3:6)-x(1:4).^2.*v(3:6).*U(1, 1:4))./(x(3:6).^3-x(1:4).^3);
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
    
    % Approxiate sf at the next step.
    % Solve sf by Sommefeld condition
    W(2, :) = W(1, :)+dt*v.*V(1, :);
    dbdyWold(2) = (24/17*W(1, end)-59/34*W(1, end-1)+4/17*W(1, end-2)+3/34*W(1, end-3))/dx+W(1, end)/x(end);
    dbdyWold(1) = 0.5*(W(1, end)-W(1, end-2))/dx+W(1, end-1)/x(end-1);
    W(2, end-1:end) = W(1, end-1:end)-dt*dbdyWold;
    
    % Correct the outer boundary of phi
    dW = my5orderdudr_Gust(W(2, :), dx);
    U(2, end-2:end) = dW(end-2:end);
    
% The second sub-step of "sneak" to nect step
elseif ICN == 2
    % Compute the next step differentiation
    dU = my5orderdudr_Gust(x.^2.*v.*U(2, :), dx);
    dU = (1./x.^2).*dU;
    dU(1) = 3*(x(2)^2*v(2)*U(2, 2)-x(1)^2*v(1)*U(2, 1))/(x(2)^3-x(1)^3);
    dU(2:5) = 3*(x(3:6).^2.*v(3:6).*U(2, 3:6)-x(1:4).^2.*v(3:6).*U(2, 1:4))./(x(3:6).^3-x(1:4).^3);
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
    W(2, end) = W(1, end)-0.5*dt*(dbdyWold(end)+(24/17*W(2, end)-59/34*W(2, end-1)+4/17*W(2, end-2)+3/34*W(2, end-3))/dx+W(2, end)/x(end));
    W(2, end-1) = W(1, end-1)-0.5*dt*(dbdyWold(end-1)+0.5*(W(2, end)-W(2, end-2))/dx+W(2, end-1)/x(end-1));
    
    % Correct the outer boundary of phi
    dW = my5orderdudr_Gust(W(2, :), dx);
    U(2, end-2:end) = dW(end-2:end);
    
    
    %----------------------------------------------------------------------
elseif ICN == 3
    dU = my5orderdudr_Gust(x.^2.*v.*U(2, :), dx);
    dU = (1./x.^2).*dU;
    dU(1) = 3*(x(2)^2*v(2)*U(2, 2)-x(1)^2*v(1)*U(2, 1))/(x(2)^3-x(1)^3);
    dU(2:5) = 3*(x(3:6).^2.*v(3:6).*U(2, 3:6)-x(1:4).^2.*v(3:6).*U(2, 1:4))./(x(3:6).^3-x(1:4).^3);
    dV = my5orderdudr_Gust(v.*V(2, :), dx);
    
    V(2, :) = V(1, :)+0.5*dt*(dU+dUold);
    U(2, :) = U(1, :)+0.5*dt*(dV+dVold);
    temp = (SOLNEARZERO*V(2, 2:3)');
    V(2, 1) = temp(2);
    U(2, 1) = 0;
    W(2, :) = W(1, :)+0.5*dt*v.*(V(1, :)+V(2, :));

    W(2, end) = W(1, end)-0.5*dt*(dbdyWold(end)+(24/17*W(2, end)-59/34*W(2, end-1)+4/17*W(2, end-2)+3/34*W(2, end-3))/dx+W(2, end)/x(end));
    W(2, end-1) = W(1, end-1)-0.5*dt*(dbdyWold(end-1)+0.5*(W(2, end)-W(2, end-2))/dx+W(2, end-1)/x(end-1));
    dW = my5orderdudr_Gust(W(2, :), dx);
    U(2, end-2:end) = dW(end-2:end);
    
    % clear
    U(1, :) = U(2, :);
    V(1, :) = V(2, :);
    W(1, :) = W(2, :);
end