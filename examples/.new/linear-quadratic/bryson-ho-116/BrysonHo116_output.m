%--------------------------------------------------------------------------
% BrysonHo116_output.m
% Output function for BrysonHo116 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo116_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% solution on T
sol(1).T = T;
sol(1).U(:,1) = BrysonHo116_U(T,in.tf,auxdata.v0,auxdata.x0);
sol(1).U(:,2) = abs(sol(1).U(:,1));
sol(1).Y = BrysonHo116_Y(T,in.tf,auxdata.v0,auxdata.x0);
sol(1).F = BrysonHo116_F(in.tf,auxdata.v0,auxdata.x0);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e5)';
    sol(2).U(:,1) = BrysonHo116_U(sol(2).T,in.tf,auxdata.v0,auxdata.x0);
    sol(2).U(:,2) = abs(sol(2).U(:,1));
    sol(2).Y = BrysonHo116_Y(sol(2).T,in.tf,auxdata.v0,auxdata.x0);
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