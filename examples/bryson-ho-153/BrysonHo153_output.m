%--------------------------------------------------------------------------
% BrysonHo153_output.m
% Output function for BrysonHo153 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo153_output(T,U,Y,P,F,p,opts)

% solution on T
sol(1).T = T;
sol(1).U = BrysonHo153_U(p.c,p.t0,p.tf,p.x0)*ones(size(T));
sol(1).Y = BrysonHo153_Y(p.c,T,p.t0,p.tf,p.x0);
sol(1).F = BrysonHo153_F(p.c,p.t0,p.tf,p.x0);

% solution on high resolution T
if opts.plotflag
    sol(2).T = linspace(p.t0,p.tf,1e4)';
    sol(2).U = BrysonHo153_U(p.c,p.t0,p.tf,p.x0)*ones(size(sol(2).T));
    sol(2).Y = BrysonHo153_Y(p.c,sol(2).T,p.t0,p.tf,p.x0);
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'Ymax';
O(2).value = max(errorU(:,1));
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(opts.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(opts.QPsolvetime);
O(5).label = 'QPsolvetime';