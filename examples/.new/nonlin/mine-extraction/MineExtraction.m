%--------------------------------------------------------------------------
% MineExtraction.m
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
obj_approach = 'string'; % 'function' or 'string' (see below)
tf = 2; % contract length
auxdata.a = 1; % profit rate
auxdata.x0 = 10; % initial ore available
auxdata.t0 = 0;
auxdata.tf = tf;
setup.auxdata = auxdata;

% time horizon
setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

% number of controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 1;

% Lagrange term
switch obj_approach
    %----------------------------------------------------------------------
    case 'string'
        setup.nonlin.lagrange = 'u1^2/y1 - a*u1';
    %----------------------------------------------------------------------
    case 'function'
        % provide function, rather than a string, for the objective function
        % NOTE: this feature is currently undocumented and undeveloped
        % 2024/04/11: currently does not work with new struct format
        setup.nonlin.lagrange = []; % only needs to have the field to work
        obj.f = { @(t,param,UYP) UYP(:,1).^2./UYP(:,2) - param(:,1).*UYP(:,1) };
        setup.internalinfo.obj = obj;
    %----------------------------------------------------------------------
end

% system dynamics
setup.lq.dynamics.A = 0;
setup.lq.dynamics.B = -1;

% symbolic data for nonlin
setup.nonlin.data.symbols = 'a';
setup.nonlin.data.values = [auxdata.a];

% simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0;% initial states
LB(1).right = 4; LB(1).matrix = auxdata.x0;
UB(2).right = 2; UB(2).matrix = auxdata.x0;% state bounds
LB(2).right = 2; LB(2).matrix = 0;
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[auxdata.x0];[auxdata.x0]];
U0 = [0;0];
setup.method.guess.X = [U0,Y0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = MineExtraction_output(T,U,Y,P,F,in,opts);

% plots
MineExtraction_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
    case 1
    % default parameters
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10; % number of nodes
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivatives = 'complex';
end

end