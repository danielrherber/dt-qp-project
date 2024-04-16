%--------------------------------------------------------------------------
% BrysonHo248_output.m
% Output function for BrysonHo248 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo248_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

[t1,t2] = BrysonHo248_T(auxdata.alpha,auxdata.beta,auxdata.gamma,in.tf);

% solution on T
sol(1).T = T;
sol(1).U = BrysonHo248_U(auxdata.alpha,auxdata.beta,auxdata.gamma,T,t1,t2,in.tf);
sol(1).Y = BrysonHo248_Y(auxdata.alpha,auxdata.beta,auxdata.gamma,T,t1,t2,in.tf);
sol(1).F = BrysonHo248_F(auxdata.alpha,auxdata.beta,auxdata.gamma,t1,t2,in.tf);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = BrysonHo248_U(auxdata.alpha,auxdata.beta,auxdata.gamma,sol(2).T,t1,t2,in.tf);
    sol(2).Y = BrysonHo248_Y(auxdata.alpha,auxdata.beta,auxdata.gamma,sol(2).T,t1,t2,in.tf);
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'Y1max';
O(2).value = max(errorY(:,2));
O(2).label = 'Y2max';
O(3).value = max(errorU(:,1));
O(3).label = 'Umax';
O(4).value = max(errorF);
O(4).label = 'F';
O(5).value = max(in.QPcreatetime);
O(5).label = 'QPcreatetime';
O(6).value = max(in.QPsolvetime);
O(6).label = 'QPsolvetime';

end