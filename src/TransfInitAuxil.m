%% Test of scalar field equatiion
%  Construct phi and pi, which is the derivative of scalar field 
%  w.r.t. space and time
% 
%  Syntax
%
%  Descriptions
%
%%
function [phicap, picap] = TransfInitAuxil_Dense(u0, drDense)
phicap = my5orderdudr_Gust(u0, drDense);
picap = zeros(1, numel(u0));
end