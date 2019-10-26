%--------------------------------------------------------------------------
% MinimumEnergyTransfer_output.m
% Output function for MinimumEnergyTransfer example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = MinimumEnergyTransfer_output(T,U,Y,P,F,in,opts)

% solution on T
[Uactual,Yactual,Factual] = MinimumEnergyTransfer_solution(T,in);
sol(1).T = T;
sol(1).U = Uactual;
sol(1).Y = Yactual;
sol(1).F = Factual;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    [Uactual,Yactual,~] = MinimumEnergyTransfer_solution(sol(2).T,in);
    sol(2).U = Uactual;
    sol(2).Y = Yactual;
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:));
O(1).label = 'Ymax';
O(2).value = max(errorU(:));
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';

end