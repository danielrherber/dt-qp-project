%--------------------------------------------------------------------------
% LQRstandard.m
% pp. 148-152 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
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
[O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts);


%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% auxiliary data
auxdata.ns = 20; % number of states
auxdata.nu = 10; % number of controls
auxdata.x0 = linspace(-5,5,auxdata.ns)'; % initial states
rng(393872382,'twister') % specific random seed
auxdata.A = sprand(auxdata.ns,auxdata.ns,0.5,1);
auxdata.B = sprand(auxdata.ns,auxdata.nu,1,1);
auxdata.R = eye(auxdata.nu);
auxdata.Q = sprand(auxdata.ns,auxdata.ns,0.2);
auxdata.Q = ((auxdata.Q)*((auxdata.Q)'))/100;
auxdata.M = 10*eye(auxdata.ns); % objective
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 10;

% system dynamics
setup.lq.dynamics.A = auxdata.A; 
setup.lq.dynamics.B = auxdata.B;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = auxdata.R; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = auxdata.Q; % x'*Q*x
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = auxdata.M; %xf'*M*xf
setup.lq.mayer = M;

% simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0; % initial states
LB(1).right = 4; LB(1).matrix = auxdata.x0; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = LQRstandard_output(T,U,Y,P,F,in,opts,setup);

% plots
LQRstandard_plot(T,U,Y,P,F,in,opts,sol)

end


%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 2;

switch num
case 1
    % default parameters
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = 0;
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 25; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
case 2
    opts.dt.defects = 'FO';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 3
    opts.dt.defects = 'PS-MI';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10; % number of intervals
    opts.dt.nn = 8; % polynomial order in each interval
end

end