%--------------------------------------------------------------------------
% Biegler9p5_output.m
% Output function for Biegler9p5 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = Biegler9p5_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% solution on T
sol(1).T = T;
sol(1).P = Biegler9p5_P();
sol(1).Y = Biegler9p5_Y(T);
sol(1).F = Biegler9p5_F();

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e5)';
    sol(2).P = Biegler9p5_P();
    sol(2).Y = Biegler9p5_Y(sol(2).T);
    sol(2).F = Biegler9p5_F();
end

% errors
errorP = abs(P-sol(1).P);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'X1max';
O(2).value = max(errorY(:,2));
O(2).label = 'X2max';
O(3).value = max(errorP);
O(3).label = 'Pmax';
O(4).value = max(errorF);
O(4).label = 'F';
O(5).value = max(in.QPcreatetime);
O(5).label = 'QPcreatetime';
O(6).value = max(in.QPsolvetime);
O(6).label = 'QPsolvetime';

end