%% SommerfeldOuterBdy
%  Correct the outer boundary of scalar field by Sommerfeld condition such
%  that there is no incoming field.
% Please refer to Arnowitt et al. (1983) and Gundlach (2004).
%
%  Syntax
%
%  Descriptions
%
%%
function Woutput = SommerfeldOuterBdy(x, W, dt, dx)
global ICN 
persistent dbdyWold

if ICN == 1
    dbdyWold(2) = (24/17*W(1, end)-59/34*W(1, end-1)+4/17*W(1, end-2)+3/34*W(1, end-3))/dx+W(1, end)/x(end);
    dbdyWold(1) = 0.5*(W(1, end)-W(1, end-2))/dx+W(1, end-1)/x(end-1);
    Woutput = W(1, end-1:end)-dt*dbdyWold;
else
    Woutput(2) = W(1, end)-0.5*dt*(dbdyWold(end)+(24/17*W(2, end)-59/34*W(2, end-1)+4/17*W(2, end-2)+3/34*W(2, end-3))/dx+W(2, end)/x(end));
    Woutput(1) = W(1, end-1)-0.5*dt*(dbdyWold(end-1)+0.5*(W(2, end)-W(2, end-2))/dx+W(2, end-1)/x(end-1));
end