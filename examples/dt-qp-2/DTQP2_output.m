%--------------------------------------------------------------------------
% DTQP2_output.m
% Output function for DTQP2 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = DTQP2_output(T,U,Y,P,F,in,opts)

% extract parameter structure
p = in.p;

% solution on T
sol(1).T = T;
sol(1).U = DTQP2_U(p.a,p.b,p.m,p.r,T,in.tf,p.omega,p.x0);
sol(1).Y = DTQP2_Y(p.a,p.b,p.m,p.r,T,in.tf,p.omega,p.x0);
sol(1).F = integral(@(t) DTQP2_L(p.a,p.b,p.m,p.r,t,in.tf,p.omega,p.x0),0,in.tf,...
    'AbsTol',eps) + p.m*sol(1).Y(end)^2;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = DTQP2_U(p.a,p.b,p.m,p.r,sol(2).T,in.tf,p.omega,p.x0);
    sol(2).Y = DTQP2_Y(p.a,p.b,p.m,p.r,sol(2).T,in.tf,p.omega,p.x0);
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