%--------------------------------------------------------------------------
% BrysonDenham_output.m
% Output function for BrysonDenham example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonDenham_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% solution on T
sol(1).T = T;
sol(1).U = BrysonDenham_U(T,auxdata.ell);
sol(1).Y = BrysonDenham_Y(T,auxdata.ell);
sol(1).F = BrysonDenham_F(auxdata.ell);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e5)';
    sol(2).U = BrysonDenham_U(sol(2).T,auxdata.ell);
    sol(2).Y = BrysonDenham_Y(sol(2).T,auxdata.ell);
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'Xmax';
O(2).value = max(errorY(:,2));
O(2).label = 'Vmax';
O(3).value = max(errorU(:,1));
O(3).label = 'Umax';
O(4).value = max(errorF);
O(4).label = 'F';
O(5).value = max(in.QPcreatetime);
O(5).label = 'QPcreatetime';
O(6).value = max(in.QPsolvetime);
O(6).label = 'QPsolvetime';

end