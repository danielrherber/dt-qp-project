%--------------------------------------------------------------------------
% BrysonHo156_output.m
% Output function for BrysonHo156 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo156_output(T,U,Y,P,F,in,opts)

% extract parameter structure
p = in.p;

% solution on T
sol(1).T = T;
sol(1).U = BrysonHo156_U(p.c,T,in.t0,in.tf,p.v0,p.omega,p.x0);
sol(1).Y = BrysonHo156_Y(p.c,T,in.t0,in.tf,p.v0,p.omega,p.x0);
sol(1).F = BrysonHo156_F(p.c,in.t0,in.tf,p.v0,p.omega,p.x0);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = BrysonHo156_U(p.c,sol(2).T,in.t0,in.tf,p.v0,p.omega,p.x0);
    sol(2).Y = BrysonHo156_Y(p.c,sol(2).T,in.t0,in.tf,p.v0,p.omega,p.x0);
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