%--------------------------------------------------------------------------
% LQRInhomogeneous.m
% pp. 175-176 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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

% tunable parameters
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
auxdata.d = LQRInhomogeneous_d(auxdata.ns); % time-varying disturbances
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 10;

% system dynamics
setup.lq.dynamics.A = auxdata.A; 
setup.lq.dynamics.B = auxdata.B;
setup.lq.dynamics.Bv = auxdata.d;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = auxdata.R/2; % u'*R*u
L(2).left = 2; L(2).right = 2; L(2).matrix = auxdata.Q/2; % x'*Q*x
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = auxdata.M/2; %xf'*M*xf
setup.lq.mayer = M;

% simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0; % initial states
LB(1).right = 4; LB(1).matrix = auxdata.x0;
setup.lq.ub = UB; setup.lq.lb = LB;


end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = LQRInhomogeneous_output(T,U,Y,P,F,in,setup,opts);

% plots
LQRInhomogeneous_plot(T,U,Y,P,F,in,opts,sol)

end


%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 1;

switch num
case 1
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 4
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-1;
end

end