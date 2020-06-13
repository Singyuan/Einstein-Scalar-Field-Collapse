function dadr = dadrODE45(a, rip, phicap, picap, rDense)
% rip : matlab compute r by tol
% r : given r in para
phicap = interp1(rDense, phicap, rip);
picap = interp1(rDense, picap, rip);

if rip == 0
    dadr = 0;
else
    dadr = (a-a^3)/(2*rip)+2*pi*a*rip*(picap^2+phicap^2);
end
end