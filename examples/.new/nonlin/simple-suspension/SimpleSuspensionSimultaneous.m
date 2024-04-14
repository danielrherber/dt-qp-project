%--------------------------------------------------------------------------
% SimpleSuspensionSimultaneous.m
% D. R. Herber and A. K. Sundarrajan, "On the uses of linear-quadratic
% methods in solving nonlinear dynamic optimization problems with direct
% transcription", in ASME International Mechanical Engineering Congress &
% Exposition, 2020.
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

% auxiliary data
auxdata = SimpleSuspensionProblemParameters;
auxdata.t0 = 0; auxdata.tf = 3;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 4;
setup.counts.np = 2;

% objective function
strL{1} = 'w1*y1^2 + ';
strL{2} = 'w2/ms^2*( p1*y2 - p2*y3 - p1*y4 + u1 )^2 + ';
strL{3} = 'w3*u1^2';
setup.nonlin.lagrange = horzcat(strL{:});

% system dynamics
str{1} =  '[';
str{end+1} = 'y2 - z0d;';
str{end+1} = '(- kt*y1 - (p1 + ct)*y2 + p2*y3 + p1*y4 + ct*z0d - u1)/mus;'; % missing y2 term
str{end+1} = 'y4 - y2;';
str{end+1} = '( p1*y2 - p2*y3 - p1*y4 + u1)/ms';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 'ct kt mus ms w1 w2 w3 z0d';
setup.nonlin.data.values = {auxdata.bt auxdata.kt auxdata.mu auxdata.ms auxdata.w1 auxdata.w2 auxdata.w3 auxdata.z0dot};

% initial state values
LB(1).right = 4;
LB(1).matrix = zeros(4,1);
UB(1).right = 4;
UB(1).matrix = zeros(4,1);

% simple parameter bounds
LB(2).right = 3;
LB(2).matrix = [auxdata.bmin;auxdata.kmin];
UB(2).right = 3;
UB(2).matrix = [auxdata.bmax;auxdata.kmax];

% rattlespace constraints
LB(3).right = 2; % states
LB(3).matrix = [-inf,-inf,-auxdata.rmax,-inf];
UB(3).right = 2; % states
UB(3).matrix = [inf,inf,auxdata.rmax,inf];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[0,0,0,0];[0,0,0,0]];
U0 = [0;0];
P0 = [LB(2).matrix';LB(2).matrix'];
setup.method.guess.X = [U0,Y0,P0];

% scaling (with knowledge of the solution)
scaling(1).right = 1; % controls
scaling(1).matrix = 400;
scaling(2).right = 2; % states
scaling(2).matrix = [0.004 0.4 auxdata.rmax 0.04];
scaling(3).right = 3; % parameters
scaling(3).matrix = [1e2;1e4];
setup.method.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = SimpleSuspension_output(T,U,Y,P,F,in,opts);

% plots
SimpleSuspension_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 1000;
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.solver.tolerance = 1e-14;
    opts.solver.maxiters = 10000;
    opts.solver.function = 'IPFMINCON';

    % opts.method.derivatives = 'internal';
    % opts.method.derivatives = 'real-forward';
    % opts.method.derivatives = 'real-central';
    opts.method.derivatives = 'complex';
    % opts.method.derivatives = 'symbolic';
end

end