%--------------------------------------------------------------------------
% DTQP_indexcolumns.m
% Create optimization variable (column) locations
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function J = DTQP_indexcolumns(a,b,c)

J = a*c+1:a*(b+c);
J = reshape(J,a,b);
J(end,:) = [];
J = J(:);

end