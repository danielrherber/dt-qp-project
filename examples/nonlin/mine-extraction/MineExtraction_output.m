%--------------------------------------------------------------------------
% MineExtraction_output.m
% Output function for Mine Extraction example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = MineExtraction_output(T,U,Y,P,F,in,opts)

% extract parameter structure
p = in.p;
a = p.a;
Tf = p.tf;
x0 = p.x0;

% solution on T
sol(1).T = T;
sol(1).U = MineExtraction_U(T,a,Tf,x0);
sol(1).Y = MineExtraction_Y(T,a,Tf,x0);
sol(1).F = MineExtraction_F(T,a,Tf,x0);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = MineExtraction_U(sol(2).T,a,Tf,x0);
    sol(2).Y  = MineExtraction_Y(sol(2).T,a,Tf,x0);
    sol(2).F = sol(1).F;
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