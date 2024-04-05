%--------------------------------------------------------------------------
% Brachistochrone_output.m
% Output function for Brachistochrone example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [O,sol] = Brachistochrone_output(T,U,Y,P,F,in,opts)

% extract parameter structure
auxdata = in.auxdata;

% failure
if isempty(P)
    O(1).value = inf;
    O(1).label = 'Umax';
    O(2).value = inf;
    O(2).label = 'Ymax';
    O(3).value = inf;
    O(3).label = 'F';
    O(4).value = max(in.QPcreatetime);
    O(4).label = 'QPcreatetime';
    O(5).value = max(in.QPsolvetime);
    O(5).label = 'QPsolvetime';
    sol = [];
    return
end

switch auxdata.variant
%----------------------------------------------------------------------
case {1,2}
% solution on T
sol(1).T = T*P(1);
[F2,Y2,U2] = Brachistochrone_solution(sol(1).T,auxdata.variant,auxdata.g,auxdata.xf,auxdata.yf);
sol(1).U = U2;
sol(1).Y = Y2;
sol(1).F = F2;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(0,F2,1e4)';
    [~,Y2,U2] = Brachistochrone_solution(sol(2).T,auxdata.variant,auxdata.g,auxdata.xf,auxdata.yf);
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

%----------------------------------------------------------------------
case 3
% solution on T
sol(1).T = T*P(1);
[F2,Y2,U2] = Brachistochrone_solution(sol(1).T,auxdata.variant,auxdata.g,auxdata.xf);
sol(1).U = U2;
sol(1).Y = Y2;
sol(1).F = F2;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(0,F2,1e4)';
    [~,Y2,U2] = Brachistochrone_solution(sol(2).T,auxdata.variant,auxdata.g,auxdata.xf);
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

%----------------------------------------------------------------------
case 4
% solution on T
sol(1).T = T*P(1);
[F2,Y2,U2] = Brachistochrone_solution(sol(1).T,auxdata.variant,auxdata.g,auxdata.xf,auxdata.theta,auxdata.h);
sol(1).U = U2;
sol(1).Y = Y2;
sol(1).F = F2;

% solution on high resolution T
if opts.general.plotflag
    sol(2).T = linspace(0,F2,1e4)';
    [~,Y2,U2] = Brachistochrone_solution(sol(2).T,auxdata.variant,auxdata.g,auxdata.xf,auxdata.theta,auxdata.h);
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