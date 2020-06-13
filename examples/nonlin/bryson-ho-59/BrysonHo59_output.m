%--------------------------------------------------------------------------
% BrysonHo59_output.m
% Output function for BrysonHo59 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = BrysonHo59_output(T,U,Y,P,F,in,opts)

% extract parameter structure
p = in.p;

% transform F
F = -F;

% solution on T
sol(1).T = T;
[U2,Y2,F2] = BrysonHo59_solution(p.tf,p.a,p.h,T);
sol(1).U = U2;
sol(1).Y = Y2;
sol(1).F = F2;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    [U2,Y2,F2] = BrysonHo59_solution(p.tf,p.a,p.h,sol(2).T);
    sol(2).U = U2;
    sol(2).Y = Y2;
    sol(2).F = F2;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorU,[],'all');
O(1).label = 'Umax';
O(2).value = max(errorY,[],'all');
O(2).label = 'Ymax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';

end