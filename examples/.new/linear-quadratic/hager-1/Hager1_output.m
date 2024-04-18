%--------------------------------------------------------------------------
% Hager1_output.m
% Output function for Hager1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = Hager1_output(T,U,Y,P,F,in,opts)

% solution on T
sol(1).T = T;
sol(1).U = Hager1_U(T);
sol(1).Y = Hager1_Y(T);
sol(1).F = Hager1_F;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = Hager1_U(sol(2).T);
    sol(2).Y = Hager1_Y(sol(2).T);
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'Xmax';
O(2).value = max(errorU(:,1));
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';

end