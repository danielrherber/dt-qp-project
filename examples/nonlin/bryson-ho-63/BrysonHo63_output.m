%--------------------------------------------------------------------------
% BrysonHo63_output.m
% Output function for BrysonHo63 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo63_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% transform F
F = -F;

% solution on T
sol(1).T = T;
[F2,E2] = BrysonHo63_solution(auxdata.tf,auxdata.V,auxdata.w,T);
sol(1).F = F2;
sol(1).E = E2;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    [F2,E2] = BrysonHo63_solution(auxdata.tf,auxdata.V,auxdata.w,sol(2).T);
    sol(2).F = F2;
    sol(2).E = E2;
end

% errors
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorF);
O(1).label = 'F';
O(2).value = max(in.QPcreatetime);
O(2).label = 'QPcreatetime';
O(3).value = max(in.QPsolvetime);
O(3).label = 'QPsolvetime';

end