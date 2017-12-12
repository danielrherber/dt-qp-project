%--------------------------------------------------------------------------
% BrysonHo166_output.m
% Output function for BrysonHo166 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo166_output(T,U,Y,P,F,p,opts)

% solution on T
sol(1).T = T;
sol(1).U = BrysonHo166_U(T,p.tf,p.v0,p.x0);
sol(1).Y = BrysonHo166_Y(T,p.tf,p.v0,p.x0);
sol(1).F = BrysonHo166_F(p.tf,p.v0,p.x0);

% solution on high resolution T
if opts.plotflag
    sol(2).T = linspace(p.t0,p.tf,1e4)';
    sol(2).U = BrysonHo166_U(sol(2).T,p.tf,p.v0,p.x0);
    sol(2).Y = BrysonHo166_Y(sol(2).T,p.tf,p.v0,p.x0);
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
O(5).value = max(opts.QPcreatetime);
O(5).label = 'QPcreatetime';
O(6).value = max(opts.QPsolvetime);
O(6).label = 'QPsolvetime';