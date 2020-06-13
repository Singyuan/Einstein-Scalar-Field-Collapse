function dalphadr = dalphadrODE45rev(alpha, rip, a, daadr, r)
% rip : matlab compute r by tol
% r : given r in para
rip = r(end)-rip;
aip = interp1(r, a, rip);
dadra = interp1(r, daadr, rip);


if rip < 1e-7
    dalphadr = 0;
else
    dalphadr = -alpha*dadra-(alpha*(aip^2-1))/rip;
end
end