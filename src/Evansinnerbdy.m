%% Evansinnerbdy
%  Correct the inner boundary of dpidr.
%  Please refer to Evans (1984).
%
%  Syntax
%
%  Descriptions
%
%%
function dU = Evansinnerbdy(x, v, U)
dU(1) = 3*(x(2)^2*v(2)*U(1, 2)-x(1)^2*v(1)*U(1, 1))/(x(2)^3-x(1)^3);
dU(2:5) = 3*(x(3:6).^2.*v(3:6).*U(1, 3:6)-x(1:4).^2.*v(3:6).*U(1, 1:4))./(x(3:6).^3-x(1:4).^3);
end