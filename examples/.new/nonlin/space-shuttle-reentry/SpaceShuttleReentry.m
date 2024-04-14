%--------------------------------------------------------------------------
% SpaceShuttleReentry.m
% pp. 247-256 of J. T. Betts, "Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming*." Society for Industrial and
% Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clearvars -except externalInput

% set options and standardize
if ~exist('externalInput','var')
    opts = localOpts;
end
DTQP_standardizedinputs2

% create setup structure
setup = createSetup;

% solve with DTQP
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% post-processing
[O,sol] = postProcessing(T,U,Y,P,F,in,opts);

%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% select test case (see below)
testnum = 2;

% test cases
switch testnum
    %----------------------------------------------------------------------
    case 1 % maximum downrange (traj09 in prbook)
    M(1).right = 5; M(1).left = 0; M(1).matrix = [0,-1,0,0,0,0]; % final states
    bankmin = deg2rad(-90); % minimum bank angle
    lonmax = deg2rad(inf); % maximum longitudinal angle
    %----------------------------------------------------------------------
    case 2 % maximum crossrange (traj21 in prbook)
    M(1).right = 5; M(1).left = 0; M(1).matrix = [0,0,-1,0,0,0]; % final states
    bankmin = deg2rad(-90); % minimum bank angle
    lonmax = deg2rad(180); % maximum longitudinal angle
    %----------------------------------------------------------------------
    case 3 % maximum crossrange with control bound (traj22 in prbook)
    M(1).right = 5; M(1).left = 0; M(1).matrix = [0,0,-1,0,0,0]; % final states
    bankmin = deg2rad(-70); % minimum bank angle
    lonmax = deg2rad(180); % maximum longitudinal angle
    %----------------------------------------------------------------------
    case 4 % maximum crossrange with heat limit (traj36 in prbook)
    % Mayer term
    M(1).right = 5; M(1).left = 0; M(1).matrix = [0,0,-1,0,0,0]; % final states
    bankmin = deg2rad(-90); % minimum bank angle
    lonmax = deg2rad(180); % maximum longitudinal angle

    % inequality constraints
    setup.nonlin.inequality.func = '[-17700*((1370823668977541*exp((7387794664292129*Re)/53592835565708902400 - (7387794664292129*(rads*y1))/53592835565708902400))/576460752303423488)^(1/2)*((7387794664292129*(vs*y4))/22517998136852480000)^(307/100)*((108842792728842046125*u1^3)/(18446744073709551616*pi^3) - (1987855047651470175*u1^2)/(288230376151711744*pi^2) + (249209700126165765*u1)/(72057594037927936*pi) - 4806323037483049/4503599627370496) - qu]';
    setup.nonlin.inequality.pathboundary = [1];
    %----------------------------------------------------------------------
end

% data
Re = 6371203.92;
cd0 = 0.07854;
cd1 = -0.35289616517697663074;
cd2 = 2.039962128348097688027;
rho0 = 1.225570827014494;
H = 7254.24;
cl0 = -0.20704;
cl1 = 1.6755577760805793917210;
S = 249.9091776;
mass = 92079.2525560557;
xmu = 3.9860319540930508801e14;
qu = 70;

% scaling
% rads = Re + 79248; vs = 7802.88;
rads = 1; vs = 1;
auxdata.rads = rads; auxdata.vs = vs; auxdata.Re = Re; % store for later

% initial conditions
h0 = 79248;
rad0 = (h0 + Re)/rads;
lon0 = deg2rad(0);
lat0 = deg2rad(0);
v0 = 7802.88/vs;
fpa0 = deg2rad(-1);
azi0 = deg2rad(90);

% final conditions
hf = 24384;
radf = (hf + Re)/rads;
vf = 762/vs;
fpaf = deg2rad(-5);

% bounds
radmin = Re/rads; radmax = rad0;
lonmin = deg2rad(-180);
latmin = deg2rad(-70); latmax = deg2rad(70);
vmin = 10/vs; vmax = 45000/vs;
fpamin = deg2rad(-80); fpamax = deg2rad(80);
azimin = deg2rad(-180); azimax = deg2rad(180);
aoamin = deg2rad(-90); aoamax = deg2rad(90);
bankmax = deg2rad(1);

% auxiliary data
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 6;
setup.counts.np = 1;

% Mayer term
setup.lq.mayer = M;

% system dynamics (see SpaceShuttleReentry_helper.m)
str{1} = '[';
str{end+1} = 'p1*((vs*y4)*sin(y5)/rads);';
str{end+1} = 'p1*(((vs*y4)*cos(y5)*sin(y6))/((rads*y1)*cos(y3)));';
str{end+1} = 'p1*(((vs*y4)*cos(y5)*cos(y6))/(rads*y1));';
str{end+1} = 'p1*((- (xmu*sin(y5))/(rads*y1)^2 - (S*rho0*(vs*y4)^2*exp((Re - (rads*y1))/H)*(cd2*u1^2 + cd1*u1 + cd0))/(2*mass))/vs);';
str{end+1} = 'p1*((cos(y5)*((vs*y4)^2/(rads*y1) - xmu/(rads*y1)^2) + (S*rho0*(vs*y4)^2*exp((Re - (rads*y1))/H)*cos(u2)*(cl0 + cl1*u1))/(2*mass))/(vs*y4));';
str{end+1} = 'p1*((((vs*y4)^2*cos(y5)*sin(y6)*tan(y3))/(rads*y1) + (S*rho0*(vs*y4)^2*exp((Re - (rads*y1))/H)*sin(u2)*(cl0 + cl1*u1))/(2*mass*cos(y5)))/(vs*y4))';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 'Re cd0 cd1 cd2 rho0 H cl0 cl1 S mass xmu rads vs qu';
setup.nonlin.data.values = [Re cd0 cd1 cd2 rho0 H cl0 cl1 S mass xmu rads vs qu];

% simple bounds
UB(1).right = 4; UB(1).matrix = [rad0,lon0,lat0,v0,fpa0,azi0]; % initial states
LB(1).right = 4; LB(1).matrix = [rad0,lon0,lat0,v0,fpa0,azi0];
UB(2).right = 5; UB(2).matrix = [radf inf inf vf fpaf inf]; % final states
LB(2).right = 5; LB(2).matrix = [radf -inf -inf vf fpaf -inf];
UB(3).right = 1; UB(3).matrix = [aoamax,bankmax]; % controls
LB(3).right = 1; LB(3).matrix = [aoamin,bankmin];
UB(4).right = 2; UB(4).matrix = [radmax lonmax latmax vmax fpamax azimax]; % states
LB(4).right = 2; LB(4).matrix = [radmin lonmin latmin vmin fpamin azimin];
UB(5).right = 3; UB(5).matrix = 4000; % parameters
LB(5).right = 3; LB(5).matrix = 0;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[rad0,lon0,lat0,v0,fpa0,azi0];[radf,lon0+deg2rad(10),lat0+deg2rad(10),vf,fpaf,deg2rad(-90)]];
U0 = [[deg2rad(0),deg2rad(0)];[deg2rad(0),deg2rad(0)]];
P0 = [1000;1000];
setup.method.guess.X = [U0,Y0,P0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = [aoamax abs(bankmin)];
scaling(2).right = 2; % states
scaling(2).matrix = [(Re + 79248) abs(lonmin) latmax (7802.88) fpamax azimax];
scaling(3).right = 3; % parameters
scaling(3).matrix = 4000;
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = SpaceShuttleReentry_output(T,U,Y,P,F,in,opts);

% plots
SpaceShuttleReentry_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.maxiters = 1000;
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-5;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 20; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.maxiters = 1000;
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-5;
end

end