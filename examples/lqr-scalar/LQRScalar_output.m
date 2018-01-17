%--------------------------------------------------------------------------
% LQRScalar_output.m
% Output function for LQRScalar example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = LQRScalar_output(T,U,Y,P,F,p,opts)

% constants
c1 = LQRScalar_C1(p.a,p.b,p.q,p.r);
c2 = LQRScalar_C2(p.a,p.b,c1,p.m,p.r,p.tf);
c3 = LQRScalar_C3(c1,c2,p.r,p.t0,p.x0);

% solution on T
sol(1).T = T;
sol(1).U = LQRScalar_U(p.a,p.b,c1,c2,c3,p.r,T);
sol(1).Y = LQRScalar_Y(c1,c2,c3,p.r,T);
sol(1).F = LQRScalar_F(p.a,p.b,c1,c2,c3,p.m,p.q,p.r,p.t0,p.tf);

% solution on high resolution T
if opts.plotflag
    sol(2).T = linspace(p.t0,p.tf,1e4)';
    sol(2).U = LQRScalar_U(p.a,p.b,c1,c2,c3,p.r,sol(2).T);
    sol(2).Y = LQRScalar_Y(c1,c2,c3,p.r,sol(2).T);
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
O(4).value = max(opts.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(opts.QPsolvetime);
O(5).label = 'QPsolvetime';