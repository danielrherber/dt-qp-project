%--------------------------------------------------------------------------
% Nonlinear1D_output.m
% Output function for Nonlinear1D example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = Nonlinear1D_output(T,U,Y,P,F,in,opts)

% solution on T
sol(1).T = T;
sol(1).U = Nonlinear1D_U(T);
sol(1).Y = Nonlinear1D_Y(T);
sol(1).F = Nonlinear1D_F;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = Nonlinear1D_U(sol(2).T);
    sol(2).Y = Nonlinear1D_Y(sol(2).T);
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorU,[],'all');
O(1).label = 'Umax';
O(2).value = max(errorY,[],'all');
O(2).label = 'Ymax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';

end