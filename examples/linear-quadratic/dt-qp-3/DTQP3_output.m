%--------------------------------------------------------------------------
% DTQP3_output.m
% Output function for DTQP3 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = DTQP3_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% solution on T
sol(1).T = T;
sol(1).U = ones(size(T))*DTQP3_U(auxdata.a1,auxdata.a2,auxdata.b,auxdata.m,auxdata.r,in.tf,auxdata.w1,auxdata.w2,auxdata.x0);
sol(1).Y = DTQP3_Y(auxdata.a1,auxdata.a2,auxdata.b,auxdata.m,auxdata.r,T,in.tf,auxdata.w1,auxdata.w2,auxdata.x0);
sol(1).F = DTQP3_F(auxdata.a1,auxdata.a2,auxdata.b,auxdata.m,auxdata.r,in.tf,auxdata.w1,auxdata.w2,auxdata.x0);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = ones(size(sol(2).T))*sol(1).U(1);
    sol(2).Y = DTQP3_Y(auxdata.a1,auxdata.a2,auxdata.b,auxdata.m,auxdata.r,sol(2).T,in.tf,auxdata.w1,auxdata.w2,auxdata.x0);
    sol(2).F = sol(1).F;
end

% errors
if auxdata.ParameterFlag
    errorU = abs(P(1)-sol(1).U);
else
    errorU = abs(F-sol(1).U);
end
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