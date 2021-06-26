%--------------------------------------------------------------------------
% LinearPendulum_output.m
% Output function for LinearPendulum example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = LinearPendulum_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% solution on T
sol(1).T = T;
sol(1).U = LinearPendulum_U(auxdata.m,auxdata.k,auxdata.umax,in.tf,T);
sol(1).Y = LinearPendulum_Y(auxdata.m,auxdata.k,auxdata.umax,auxdata.x0,auxdata.v0,in.tf,T);
sol(1).F = LinearPendulum_F(auxdata.m,auxdata.k,auxdata.umax,auxdata.x0,auxdata.v0,in.tf);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = LinearPendulum_U(auxdata.m,auxdata.k,auxdata.umax,in.tf,sol(2).T);
    sol(2).Y = LinearPendulum_Y(auxdata.m,auxdata.k,auxdata.umax,auxdata.x0,auxdata.v0,in.tf,sol(2).T);
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