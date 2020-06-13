function daadr = daadrexplicit(aAdap, Oldr, rAdap, phicap, picap)
% Doesn's update OldrAdap
phicapAdap = interp1(Oldr, phicap, rAdap);
picapAdap = interp1(Oldr, picap, rAdap);

daadr = (1-aAdap.^2)./(2.*rAdap)+2*pi.*rAdap.*(phicapAdap.^2+picapAdap.^2);
daadr(1) = 0;
end