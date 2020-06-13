%% solvea
% Solve the ratio of coordinate time and proper time, i.e. lapse function 
% "alpha"
% Assume the metric is ds = alpha^2 dt^2 + a^2 dr^2 + ...
%
%  Syntax
%
%  Descriptions
%
%%
function alphaAdap = solvealpha(rAdap, aAdap, phicap, picap, odeinput4alpha, rDense)
global ICN
if ICN == 1
    % Compute the da/dr
    daadr = daadrexplicit(aAdap, rDense, rAdap, phicap, picap);
    
    % Solve alpha from outside. Hence, the initial condition of alpha is
    % end of a
    icalpha = aAdap(end);
    [rAdapsec, alphaAdap] = ode45(@(rip, alpha) dalphadrODE45rev(alpha, rip, aAdap, daadr, rAdap), odeinput4alpha.range, icalpha, odeinput4alpha.opts);
    
    % 
    alphaAdap = interp1(fliplr(odeinput4alpha.range(end)-rAdapsec'), fliplr(alphaAdap'), rAdap);

    if isnan(alphaAdap(1))
        disp('warning alpha')
        XI = interp1(rAdap(4:end), alphaAdap(4:end), rAdap(1:3), 'pchip','extrap');
        alphaAdap(1:3) = XI;
    end
else
    daadr = daadrexplicit(aAdap, rDense, rAdap, phicap, picap);
    % use ode45 to compute the alpha
    icalpha = aAdap(end);
    [rAdapsec, alphaAdap] = ode23(@(rip, alpha) dalphadrODE45rev(alpha, rip, aAdap, daadr, rAdap), odeinput4alpha.range, icalpha, odeinput4alpha.opts);
    alphaAdap = interp1(fliplr(odeinput4alpha.range(end)-rAdapsec'), fliplr(alphaAdap'), rAdap);
    if isnan(alphaAdap(1))
        disp('warning alpha')
        XI = interp1(rAdap(4:end), alphaAdap(4:end), rAdap(1:3), 'pchip','extrap');
        alphaAdap(1:3) = XI;
    end
end
end