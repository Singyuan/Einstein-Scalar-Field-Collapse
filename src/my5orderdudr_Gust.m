%% my5orderdudr_Gust
%  Compute the derivative by summation by part.
% Please refer to Gustanfsson (2008) and Strand (1994)
%
%  Syntax
%
%  Descriptions
%
%%
function dudr = my5orderdudr_Gust(u, drDense)
uDiff = u(3:end)-u(1:end-2);
dudr2h = uDiff./(2*drDense);
dudr2h = [0 dudr2h 0];

uDiff = u(5:end)-u(1:end-4);
dudr4h = uDiff./(4*drDense);
dudr4h = [0 0 dudr4h 0 0];

dudr = (4*dudr2h-1*dudr4h)/3;

dudr(1) = (-24/17*u(1)+59/34*u(2)-4/17*u(3)-3/34*u(4))/drDense;
dudr(2) = (-1/2*u(1)+1/2*u(3))/drDense;
dudr(3) = (4/43*u(1)-59/86*u(2)+59/86*u(4)-4/43*u(5))/drDense;
dudr(4) = (3/98*u(1)-59/98*u(3)+32/49*u(5)-4/49*u(6))/drDense;

dudr(end) = (24/17*u(end)-59/34*u(end-1)+4/17*u(end-2)+3/34*u(end-3))/drDense;
dudr(end-1) = (1/2*u(end)-1/2*u(end-2))/drDense;
dudr(end-2) = (-4/43*u(end)+59/86*u(end-1)-59/86*u(end-3)+4/43*u(end-4))/drDense;
dudr(end-3) = (-3/98*u(end)+59/98*u(end-2)-32/49*u(end-4)+4/49*u(end-5))/drDense;
end


