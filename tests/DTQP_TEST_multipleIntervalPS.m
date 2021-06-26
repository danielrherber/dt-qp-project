%--------------------------------------------------------------------------
% DTQP_TEST_multipleIntervalPS.m
% Test multiple-interval pseudospectral method approach that is under
% development
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:5;
% tests = 4;

% go through the tests
for k = 1:length(tests)

    clear opts

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        t0 = 0; tf = 20;
        ni = 20;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 4;
        %------------------------------------------------------------------
        case 2
        t0 = 0; tf = 20;
        ni = 4;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 7;
        %------------------------------------------------------------------
        case 3
        t0 = 0; tf = 20;
        ni = 3;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 4;
        %------------------------------------------------------------------
        case 4
        t0 = 0; tf = 20;
        ni = 40;
        xT = linspace(t0,tf,ni+1);
        opts.dt.nt = 7;
        %------------------------------------------------------------------
        case 5
        t0 = 0; tf = 20;
        ni = 4;
        xT = linspace(t0,tf,ni+1);
        opts.dt(1).nt = 4;
        opts.dt(2).nt = 5;
        opts.dt(3).nt = 6;
        opts.dt(4).nt = 10;
    end

    % run the test and time
    [setup,opts] = problem(opts,xT);
    [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

    % test analysis
    flag = 'preliminary'; closeflag = false; DTQP_plotCommon;
    hf = figure(1); hf.Color = wcolor;
    plot(T,U,'linewidth',2); hold on
    flag = 'axis'; DTQP_plotCommon;

    hf = figure(2); hf.Color = wcolor;
    plot(T,Y,'linewidth',2); hold on
    flag = 'axis'; DTQP_plotCommon;

    hf = figure(3); hf.Color = wcolor;
    plot(T,Y(:,1),'linewidth',2); hold on
    flag = 'axis'; DTQP_plotCommon;

end

% problem structure
function [setup,opts] = problem(opts,xT)

% tunable parameters
auxdata.ell = 1/9;

% options
[opts.dt.defects] = deal('PS');
[opts.dt.quadrature] = deal('G');
[opts.dt.mesh] = deal('LGL');

%% setup
% setup = BrysonHo166_local;
% setup = GreenhouseClimate_local;
setup = MinimumEnergyTransfer_local;

setup_temp = setup;
UB = setup_temp.UB;
LB = setup_temp.LB;

setup(1:(length(xT)-1)) = setup_temp;

%% linkage equality constraints
idx = 0;

% states
q = eye(length(setup_temp.A));
for k = 1:length(setup_temp.A)
    idx = idx + 1;
    LY(idx).left.linear.right = 5; % final states
    LY(idx).left.linear.matrix = q(:,k);
    LY(idx).right.linear.right = 4; % initial states
    LY(idx).right.linear.matrix = -q(:,k);
    LY(idx).b = 0;
end

% controls
q = eye(length(setup_temp.B));
for k = 1:size(setup_temp.B,2)
    idx = idx + 1;
    LY(idx).left.linear.right = 7; % final controls
    LY(idx).left.linear.matrix = q(:,k);
    LY(idx).right.linear.right = 6; % initial controls
    LY(idx).right.linear.matrix = -q(:,k);
    LY(idx).b = 0;
end

%% phases
% combine structures
for phs = 1:length(xT)-1

    % extract
    UB_right = UB;
    LB_right = LB;

    % remove 4 if not first phase
    if phs ~= 1
        UB_right([UB_right.right] == 4) = [];
        LB_right([LB_right.right] == 4) = [];
    end

    % remove 5 if not final phase
    if phs ~= length(xT)-1
        UB_right([UB_right.right] == 5) = [];
        LB_right([LB_right.right] == 5) = [];
        setup(phs).M = [];
    end

    % assign
    setup(phs).UB = UB_right;
    setup(phs).LB = LB_right;

    % mesh boundary
    setup(phs).t0 = xT(phs);
    setup(phs).tf = xT(phs+1);

    % parameters
    setup(phs).auxdata = auxdata;

    % linkage constraints
    if phs < length(xT)-1
        setup(phs).LY = LY;
    end

end

end

function setup = BrysonHo166_local

%% tunable parameters
tf = 20; % time horizon
auxdata.x0 = -0.5; auxdata.v0 = 1; % other

%% setup
t0 = 0;

% system dynamics
A = [0 1;-1 0];
B = [0;1];

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1/2; % 1/2*u^2

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = [auxdata.x0;auxdata.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [auxdata.x0;auxdata.v0];

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = [0;0];
UB(2).right = 5; % final states
UB(2).matrix = [0;0];

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

end

function setup = GreenhouseClimate_local

% problem parameters
xtf = 48;
xp1 = 7.5e-8;
xp2 = 1;
xp3 = 0.1;
xp4 = 4.55e-4;
xp5 = 136.4;
umax = 100;

% system dynamics
A = cell(2);
A{1,1} = 0;
A{1,2} = @(t) xp1*( max(0,800*sin(4*pi*t/xtf - 0.65*pi)) );
A{2,1} = 0;
A{2,2} = -xp2;
setup.A = A;
setup.B = [0;xp3];
d = cell(2,1);
d{2,1} = @(t) xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi));
setup.d = d;


auxdata.xtf = xtf; auxdata.xp1 = xp1; auxdata.xp2 = xp2; auxdata.xp3 = xp3;
auxdata.xp4 = xp4; auxdata.xp1 = xp5;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2;

% time horizon
auxdata.t0 = 0; auxdata.tf = xtf;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-xp5,0];

% Lagrange term
L(1).left = 0; L(1).right = 1; L(1).matrix = xp4;

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,10]; % initial states
LB(1).right = 4; LB(1).matrix = [0,10];
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;

% combine structures
setup.M = M; setup.L = L; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

end

function setup = MinimumEnergyTransfer_local

%% tunable parameters
t0 = 0; tf = 2; % time horizon
ny = 18; % number of states
nu = 6; % number of controls

% system dynamics
rng(83233683,'twister') % random number seed
Adensity = rand;
Aeig = -2 + (2 - -2).*rand(ny,1);
auxdata.A = sprandsym(ny,Adensity,Aeig);
auxdata.B = -10 + (10 - -10).*rand(ny,nu);

% initial states
auxdata.y0 = 10*rand(ny,1);
auxdata.yf = zeros(ny,1);

% check controllability
try
	Co = ctrb(auxdata.A,auxdata.B);
    if rank(Co) ~= ny
        warning('system is not controllable')
    end
catch
    warning('unable to check controllability')
end

%% setup
% system dynamics
A = auxdata.A;
B = auxdata.B;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = eye(nu);

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = auxdata.y0;
UB(1).right = 4; % initial states
UB(1).matrix = auxdata.y0;

% final conditions
LB(2).right = 5; % final states
LB(2).matrix = auxdata.yf;
UB(2).right = 5; % final states
UB(2).matrix = auxdata.yf;

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.auxdata = auxdata;

end