%--------------------------------------------------------------------------
% SimpleSASA_output.m
% Output function for SimpleSASA example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = SimpleSASA_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;
J = auxdata.J;
umax = auxdata.umax;
f = auxdata.tf;

% transform F
F = abs(F);

% solution on T
sol(1).T = T;
sol(1).U = SimpleSASA_U(T,umax,f);
sol(1).Y = [SimpleSASA_Y1(T,umax,f,J),SimpleSASA_Y2(T,umax,f,J)];
sol(1).P = SimpleSASA_P(f,J);
sol(1).F = SimpleSASA_F(umax,f,J);

% solution on high resolution T
if opts.general.plotflag
    Ts = SimpleSASA_Ts(f);
    sol(2).T = unique([linspace(in.t0,in.tf,1e4)';Ts-eps;Ts+eps]);
    sol(2).U = SimpleSASA_U(sol(2).T,umax,f);
    sol(2).Y  = [SimpleSASA_Y1(sol(2).T,umax,f,J),SimpleSASA_Y2(sol(2).T,umax,f,J)];
    sol(2).P = sol(1).P;
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorP = abs(P-sol(1).P);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorU,[],'all');
O(1).label = 'Umax';
O(2).value = max(errorY,[],'all');
O(2).label = 'Ymax';
O(3).value = max(errorP,[],'all');
O(3).label = 'P';
O(4).value = max(errorF);
O(4).label = 'F';
O(5).value = max(in.QPcreatetime);
O(5).label = 'QPcreatetime';
O(6).value = max(in.QPsolvetime);
O(6).label = 'QPsolvetime';

end