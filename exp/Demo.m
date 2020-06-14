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
% Global variabvle about iteration of CN
global ICN

% Set local variables
myParameters(); % set all parameter
SetVariable()
load('para.mat')
load('vari.mat')

% Initial step t = 0
totaltime = 0.0;
counter = 1;
plottime = 0.0;
recordtime = 0.0;
recorder = 1;
iter = 1;
sf = SetInitSf(rDense);
sf.origin(1, :) = InitCond(rDense, phi0, r0, pow, delta, q);
[sf.phi(1, :), sf.pi(1, :)] = TransfInitAuxil(sf.origin(1, :), drDense);

% create movie or plot at initial step
try 
    VideoMaker(rDense, sf.origin(1, :), totaltime)
    frame(counter) = getframe(gcf);
    counter = counter+1;
    plottime = plottime+0.2;
    ismovie = 1;
catch
    FigureMaker(rDense, sf.origin(1, :), totaltime, counter)
    counter = counter+1;
    plottime = plottime+1.0;
    ismovie = 0;
end


% main loop
while totaltime < T+0.2
    % The first sub-step approximate the next time step
    ICN = 1;
    % solve a
    [rAdap, aAdap] = solvea(sf.phi(1, :), sf.pi(1, :), rDense, odeinput4a);
    
    % solve alpha
    alphaAdap = solvealpha(rAdap, aAdap, sf.phi(1, :), sf.pi(1, :), odeinput4a, rDense);
    
    % create movie or plot
    if totaltime > plottime
        if ismovie == 1
            VideoMaker(rDense, sf.origin(1, :), totaltime)
            frame(counter) = getframe(gcf);
            counter = counter+1;
            plottime = plottime+0.2;
        else
            FigureMaker(rDense, sf.origin(1, :), totaltime, counter)
            counter = counter+1;
            plottime = plottime+1.0;
        end
    end
    
    % Record the data
    if totaltime >= recordtime
        data = DataPlot(data, recorder, totaltime, alphaAdap);
        recorder = recorder+1;
        recordtime = recordtime+0.2;
    end
    
    % prepare the next step Auxil coordinate
    dt = CFLcond(aAdap, alphaAdap, drDense, fg);
    v = interp1(rAdap, alphaAdap./aAdap, rDense);
    [sf.phi, sf.pi, sf.origin] = IterCrankNicolson(sf.phi, sf.pi, sf.origin, v, rDense, drDense, dt);
    totaltime = totaltime+dt;
    
    % The second sub-step to "sneak" next step
    ICN = 2;
    % solve a
    [rAdap, aAdap] = solvea(sf.phi(2, :), sf.pi(2, :), rDense, odeinput4a);
    
    % solve alpha
    alphaAdap = solvealpha(rAdap, aAdap, sf.phi(2, :), sf.pi(2, :), odeinput4a, rDense);
    
    % prepare the next step Auxil coordinate
    dt = CFLcond(aAdap, alphaAdap, drDense, fg);
    v = interp1(rAdap, alphaAdap./aAdap, rDense);
    [sf.phi, sf.pi, sf.origin] = IterCrankNicolson(sf.phi, sf.pi, sf.origin, v, rDense, drDense, dt);
    
    % The third sub-step to correct the next step
    ICN = 3;
    % solve a
    [rAdap, aAdap] = solvea(sf.phi(2, :), sf.pi(2, :), rDense, odeinput4a);
    
    % solve alpha
    alphaAdap = solvealpha(rAdap, aAdap, sf.phi(2, :), sf.pi(2, :), odeinput4a, rDense);
    
    % prepare the next step Auxil coordinate
    dt = CFLcond(aAdap, alphaAdap, drDense, fg);
    v = interp1(rAdap, alphaAdap./aAdap, rDense);
    [sf.phi, sf.pi, sf.origin] = IterCrankNicolson(sf.phi, sf.pi, sf.origin, v, rDense, drDense, dt);
    
    if mod(iter, 20) == 0
        disp(['Amplitude:' num2str(phi0), ', time: ', num2str(totaltime), ',  dt: ', num2str(dt), ',  j: ', num2str(iter)])
    end
    iter = iter+1;
end

% Output the movie or plot
if ismovie == 1
    video = VideoWriter([pwd '\..\result\scalarfield\test.avi']);
    video.FrameRate = 5;
    open(video)
    writeVideo(video, frame);
    close(video);
end

%Output the data
h = plot(data.timecoord, data.centralalpha,'LineWidth',1.5);
xlabel('time')
axis([-inf, inf, 0, 1])
title('central of alpha')
drawnow
filename = '\..\result\data\centalalpha';
saveas(h, [pwd filename], 'png');