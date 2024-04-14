%--------------------------------------------------------------------------
% SimpleSuspensionProblemParameters.m
% Problem parameters for SimpleSuspension example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function p = SimpleSuspensionProblemParameters

% suspension parameters
p.bt = 0;
p.kt = 232.5e3;
p.mu = 65;
p.ms = 325;

% rattlespace constraint
p.rmax = 0.04;

% objective weights
p.w1 = 1e5;
p.w2 = 0.5;
p.w3 = 1e-5;

% plant variable upper and lower bounds
p.bmin = 1e2; % damper min
p.bmax = 1e5; % damper max
p.kmin = 1e2; % spring min
p.kmax = 1e6; % spring max

% road input
p.v = 10;
load('IRI_737b','road_x','road_z');
road_t = road_x./p.v;
PP = spline(road_t,road_z);
PPd = fnder(PP,1);
p.z0dot = @(t) ppval(PPd,t);
p.z0 = @(t) ppval(PP,t);

end