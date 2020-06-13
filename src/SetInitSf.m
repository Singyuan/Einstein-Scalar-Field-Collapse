%% SetInitSf
%  Set initial scalar field structure include of derivative w.r.t space and
%  time
%
%  Syntax
%
%  Descriptions
%
%%
function sf = SetInitSf(rDense)
sf.phi = zeros(2, numel(rDense));
sf.pi = zeros(size(sf.phi));
sf.origin = zeros(size(sf.phi));
end