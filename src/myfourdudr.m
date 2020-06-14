%% myfourdudr
%  Construct the diffusion term by 4 times derivative.
%  Please refer to Kreiss and Oliger (1973)
%
%  Syntax
%
%  Descriptions
%
%%
function DD = myfourdudr(u)
uDiff = u(5:end)+u(1:end-4)-4*u(4:end-1)-4*u(2:end-3)+6*u(3:end-2);
dudr = uDiff;
DD = [0 0 dudr 0 0];
DD(2) = -(u(1)-2*u(2)+u(3));
DD(3) = -(u(2)-2*u(3)+u(4));
DD(4) = -(u(3)-2*u(4)+u(5));
DD(end-1) = -(u(end-2)-2*u(end-1)+u(end));
DD(end-2) = -(u(end-3)-2*u(end-2)+u(end-1));
DD(end-3) = -(u(end-4)-2*u(end-3)+u(end-2));
end


