%--------------------------------------------------------------------------
% MultiphaseParameter_output.m
% Output function for MultiphaseParameter example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = MultiphaseParameter_output(T,U,Y,P,F,in,opts)

% solution on T
[Yopt,Popt,Fopt] = MultiphaseParameter_solution(in);
sol(1).T = T;
sol(1).P = Popt;
sol(1).Y = Yopt;
sol(1).F = Fopt;

% solution on high resolution T
if opts.general.plotflag
    % initialize
    IN = in; T2 = []; T3 = [];

    % construct high resolution T
    T1 = linspace(in(1).t0,in(1).tf,1e4)';
    IN(1).t = T1;
    if length(in) > 1
        T2 = linspace(in(2).t0,in(2).tf,1e4)';
        IN(2).t = T2;
    end
    if length(in) > 2
        T3 = linspace(in(3).t0,in(3).tf,1e4)';
        IN(3).t = T3;
    end

    % obtain the solution
    [Yopt,~,~] = MultiphaseParameter_solution(IN);

    % combine
    sol(2).T = [T1;T2;T3];
    sol(2).P = sol(1).P;
    sol(2).Y = Yopt;
    sol(2).F = sol(1).F;
end

% errors
errorP = abs(P-sol(1).P);
errorY = abs(Y-sol(1).Y);
errorF = abs(F-sol(1).F);

% outputs
O(1).value = max(errorY(:));
O(1).label = 'Ymax';
O(2).value = max(errorP);
O(2).label = 'Umax';
O(3).value = max(errorF);
O(3).label = 'F';
O(4).value = max(in(1).QPcreatetime);
O(4).label = 'QPcreatetime';
O(5).value = max(in(1).QPsolvetime);
O(5).label = 'QPsolvetime';

end