%--------------------------------------------------------------------------
% NeuenhofenKerriganX2_output.m
% Output function for NeuenhofenKerriganX2 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = NeuenhofenKerriganX2_output(T,U,Y,P,F,in,opts)

% solution on T
[To,Uo,Yo,Fo] = NeuenhofenKerriganX2_solution(T);

% assign
sol(1).T = To;
sol(1).U = Uo;
sol(1).Y = Yo;
sol(1).F = Fo;

% solution on high resolution T
if opts.general.plotflag

    % solution on T
    [To,Uo,Yo,Fo] = NeuenhofenKerriganX2_solution(linspace(in.t0,in.tf,1e4)');

    % assign
    sol(2).T = To;
    sol(2).U = Uo;
    sol(2).Y = Yo;
    sol(2).F = Fo;

end

% errors
errorU = max(abs(U-sol(1).U),[],'all');
errorY = max(abs(Y-sol(1).Y),[],'all');
errorF = max(abs(F-sol(1).F),[],'all');

% outputs
O(1).value = max(errorU);
O(1).label = 'U';
O(2).value = max(errorY);
O(2).label = 'Y';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in.QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in.QPsolvetime);
O(5).label = 'QPsolvetime';

end