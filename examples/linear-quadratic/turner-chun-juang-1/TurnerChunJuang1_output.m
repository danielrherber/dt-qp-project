%--------------------------------------------------------------------------
% TurnerChunJuang1_output.m
% Output function for TurnerChunJuang1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = TurnerChunJuang1_output(T,U,Y,P,F,in,opts)

% extract parameter structure
p = in.p;

% solution on T
sol(1).T = T;
sol(1).U = TurnerChunJuang1_U(p.S,p.a,p.em,p.eta,in.tf,T,in.t0,p.um,p.y0);
sol(1).Y = TurnerChunJuang1_Y(p.S,p.a,p.em,p.eta,in.tf,T,in.t0,p.um,p.y0);
sol(1).F = TurnerChunJuang1_F(p.S,p.a,p.em,p.eta,in.tf,in.t0,p.um,p.y0);

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(in.t0,in.tf,1e4)';
    sol(2).U = TurnerChunJuang1_U(p.S,p.a,p.em,p.eta,in.tf,sol(2).T,in.t0,p.um,p.y0);
    sol(2).Y = TurnerChunJuang1_Y(p.S,p.a,p.em,p.eta,in.tf,sol(2).T,in.t0,p.um,p.y0);
    sol(2).F = sol(1).F;
end

% errors
errorU = abs(U-sol(1).U);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:,1));
O(1).label = 'Ymax';
O(2).value = max(errorU(:,1));
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';

end