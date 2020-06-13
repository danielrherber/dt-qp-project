%--------------------------------------------------------------------------
% GasAbsorber_output.m
% Output function for GasAbsorber example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = GasAbsorber_output(T,U,Y,P,F,in,opts,setup)

% solution
sol = []; % no exact solution

% outputs
O(1).value = max(in.QPcreatetime);
O(1).label = 'QPcreatetime';
O(2).value = max(in.QPsolvetime);
O(2).label = 'QPsolvetime';

end