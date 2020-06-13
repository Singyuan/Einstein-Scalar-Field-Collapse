%% Test of scalar field equatiion
%  Test of scalar field equatiion
%  Please refer to https://singyuan.github.io/projects/NTU/Numerical_Relativity.html
%
%  Syntax
%
%  Descriptions
%
%%
clear;clc;close all;
addpath('./../src')

% Set global variables
global ICN % global variabvle about iteration of CN
myParameters(); % set all parameter
load('para.mat')
sf = SetInitSf(rDense);

% Set the ode solver
odeinput4a.opts = odeset('RelTol', 1e-5, 'MaxStep', MaxStep4a, 'Events', @StopEvents);
odeinput4a.range = [0, rend];
odeinput4a.ic = ica;
odeinput4alpha.opts = odeset('RelTol', 1e-5, 'MaxStep', MaxStep4alpha);
odeinput4alpha.range = [0, rend];

% Initial step t = 0
totaltime = 0.0;
counter = 0;
iter = 1;
sf.origin(1, :) = InitCond(rDense, phi0, r0, pow, delta, q);
[sf.phi(1, :), sf.pi(1, :)] = TransfInitAuxil(sf.origin(1, :), drDense);

% 
while totaltime < T
    ICN = 1;
    % solve a
    [rAdap, aAdap] = solvea(sf.phi(1, :), sf.pi(1, :), rDense, odeinput4a);
    
    % solve alpha
    alphaAdap = solvealpha(rAdap, aAdap, sf.phi(1, :), sf.pi(1, :), odeinput4a, rDense);
    
    if mod(iter, 30) == 1
        plot([-fliplr(rDense(2:end)) rDense], [fliplr(sf.origin(1, 2:end)) sf.origin(1, :)],'LineWidth',1.5)
        xlabel('r')
        axis([-inf, inf, -0.15, 0.15])
        tit = ['time = ', num2str(totaltime)];
        title(tit)
        drawnow
        filname = [num2str(counter), '.png'];
%         frame(counter) = getframe(gcf);
        counter = counter+1;
    end
    
    % prepare the next step Auxil coordinate
    dt = CFLcond(aAdap, alphaAdap, drDense, fg);
    v = interp1(rAdap, alphaAdap./aAdap, rDense);
    [sf.phi, sf.pi, sf.origin] = IterCrankNicolson(sf.phi, sf.pi, sf.origin, v, rDense, drDense, dt);
    totaltime = totaltime+dt;
    
    ICN = 2;
    % solve a
    [rAdap, aAdap] = solvea(sf.phi(2, :), sf.pi(2, :), rDense, odeinput4a);
    
    % solve alpha
    alphaAdap = solvealpha(rAdap, aAdap, sf.phi(2, :), sf.pi(2, :), odeinput4a, rDense);
    
    % prepare the next step Auxil coordinate
    dt = CFLcond(aAdap, alphaAdap, drDense, fg);
    v = interp1(rAdap, alphaAdap./aAdap, rDense);
    [sf.phi, sf.pi, sf.origin] = IterCrankNicolson(sf.phi, sf.pi, sf.origin, v, rDense, drDense, dt);
    
    ICN = 3;
    % solve a
    [rAdap, aAdap] = solvea(sf.phi(2, :), sf.pi(2, :), rDense, odeinput4a);
    
    % solve alpha
    alphaAdap = solvealpha(rAdap, aAdap, sf.phi(2, :), sf.pi(2, :), odeinput4a, rDense);
    
    % prepare the next step Auxil coordinate
    dt = CFLcond(aAdap, alphaAdap, drDense, fg);
    v = interp1(rAdap, alphaAdap./aAdap, rDense);
    [sf.phi, sf.pi, sf.origin] = IterCrankNicolson(sf.phi, sf.pi, sf.origin, v, rDense, drDense, dt);
    
    if mod(iter, 5) == 0
        disp(['Amplitude:' num2str(phi0), ', time: ', num2str(totaltime), ',  dt: ', num2str(dt), ',  j: ', num2str(iter)])
    end
    iter = iter+1;
end