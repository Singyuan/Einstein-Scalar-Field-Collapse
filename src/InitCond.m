%% InitCond
%  Construct initial scalar field
%
%  Syntax
%
%  Descriptions
%
%%
function u0 = InitCond(rDense, phi0, r0, pow, delta, q)
u0 = phi0.*rDense.^pow.*exp(-((rDense-r0)./delta).^q);
end
