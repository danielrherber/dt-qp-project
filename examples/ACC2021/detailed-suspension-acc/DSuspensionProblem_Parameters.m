%--------------------------------------------------------------------------
% SimpleDSuspensionProblem_Parameters.m
% Problem parameters for Detailed Suspension example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function p = DSuspensionProblem_Parameters(p)

% final time
p.tf = 2;

% limits for plant design
p.xpLowerBound = [0.005 0.05 0.02 3 0.003 0.03 0.1];
p.xpUpperBound = [0.02 0.4 0.5 16 0.012 0.08 0.3];

% suspension parameters
p.bt = 0;
p.kt = 232.5e3;
p.mu = 65;
p.ms = 325;

% rattlespace constraint
p.rmax = 0.04;
p.smax = 0.04;

% objective weights
p.w1 = 1e5;
p.w2 = 0.5;
p.w3 = 1e-5;
p.wramp = 0.01;
p.wrough = 1;

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

% ramp input
p.ramp_in = 25/100*10;

% mesh from other studies
% p.Allison2014b.rough = load('RoughMeshAllison2014b','t');

% spring and damper design parameters
p.g = 9.81;               % acceleration due to gravity
p.G = 77.2e+9;            % shear modulus
p.Q = 1.75;
p.rho1 = 850;             % density of oil
p.kv = 7500;              % damper valve spring rate(N/m)
p.Afa = 0.1;              % damper valve area factor
p.eta = 0.9;              % damper valve circumference factor
p.Cd = 0.7;               % valve discharge coefficient
p.Pallow = 4.75e6;        % maximum allowed damper pressure(pa)
p.x3d_allow = 5.0;        % maximum damper velocity(m/s)
p.xv_allow = 0.03;        % maximum damper valve lift(m)
p.nd = 1.2;               % design safety factor
p.m = 0.108;
p.A = 1974;               % tensile strength
p.L0max = 0.40;           % specified pocket length
p.DoMax = 0.25;           % outer spring diameter
p.dsc = 0.009;            % clearance
p.dwt = 0.002;            % damper wall thickness
p.ld1 = 0.02;             % Space required for damper above and below piston range
p.LB = 0.02;              % Space required for damper above and below piston range
p.ld2 = 0.04;             % Space required for damper above and below piston range
p.ld3 = 0.02;             % Space required for damper above and below piston range
p.xi = 0.9;               % Limit of exposed outer circumference

end