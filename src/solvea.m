%% solvea
% Solve the ratio of coordinate length and proper length, i.e. "a"
% Assume the metric is ds = alpha^2 dt^2 + a^2 dr^2 + ...
%
%  Syntax
%
%  Descriptions
%
%%
function [rAdap, aAdap] = solvea(phicap, picap, rDense, odeinput4a)
global ICN
if ICN == 1
    [rAdap, aAdap] = ode45(@(rip, a) dadrODE45(a, rip, phicap, picap, rDense), odeinput4a.range, odeinput4a.ic, odeinput4a.opts);
    rAdap = rAdap';
    aAdap = aAdap';
    
    % Correct the initial condtion
    aAdap(1) = 1;
else
    [rAdap, aAdap] = ode23(@(rip, a) dadrODE45(a, rip, phicap, picap, rDense), odeinput4a.range, odeinput4a.ic, odeinput4a.opts);
    rAdap = rAdap';
    aAdap = aAdap';
    
    % Correct the initial condtion
    aAdap(1) = 1;
end

% Error if it cannot solve
if (rAdap(end) ~= odeinput4a.range(end))||isnan(aAdap(end))
    error('Good Game')
end
end