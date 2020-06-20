%--------------------------------------------------------------------------
% DTQP_DEFECTS_index_columns.m
% Create optimization variable (column) locations
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function J = DTQP_DEFECTS_index_columns(a,b,c)

J = a*c+1:a*(b+c);
J = reshape(J,a,b);
J(end,:) = [];
J = J(:);

end