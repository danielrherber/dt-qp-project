%--------------------------------------------------------------------------
% BettsBiehnCampbell1_output.m
% Output function for BettsBiehnCampbell1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BettsBiehnCampbell1_output(T,U,Y,P,F,in,opts)

% solution on T
sol(1).T = T;
sol(1).U = BettsBiehnCampbell1_U(T);
sol(1).Y = BettsBiehnCampbell1_Y(T);
sol(1).F = BettsBiehnCampbell1_F;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = BettsBiehnCampbell1_U(sol(2).T);
    sol(2).Y = BettsBiehnCampbell1_Y(sol(2).T);
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

end