%--------------------------------------------------------------------------
% LinearPendulum.m
% Based on the linearized pendulum example in:
% A. Bressan, "Viscosity Solutions of Hamilton-Jacobi Equations and Optimal
% Control Problems," Lecture Notes, pp. 30-31, 2011.
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
tf = 5;         auxdata.m = 1/6;  auxdata.k = 1; auxdata.x0 = -4;
auxdata.v0 = 2; auxdata.umax = 2;
setup.auxdata = auxdata;

% original parameters
% tf = 20; auxdata.m = 1; auxdata.k = 1;  auxdata.x0 = 0;
% auxdata.v0 = 0;  auxdata.umax = 1;

% time horizon
setup.t0 = 0; setup.tf = tf;

% system dynamics
setup.lq.dynamics.A = [0 1;-auxdata.k/auxdata.m 0];
setup.lq.dynamics.B = [0;1/auxdata.m];

% Mayer term
M(1).right = 5; M(1).matrix = -[1 0]; % final states
setup.lq.mayer = M;

% simple bounds
LB(1).right = 4; % initial states
LB(1).matrix = [auxdata.x0 auxdata.v0];
UB(1).right = 4; % initial states
UB(1).matrix = [auxdata.x0 auxdata.v0];
LB(2).right = 1; % controls
LB(2).matrix = -auxdata.umax;
UB(2).right = 1; % controls
UB(2).matrix = auxdata.umax;
setup.lq.ub = UB; setup.lq.lb = LB;


end



%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = LinearPendulum_output(T,U,Y,P,F,in,opts);

% plots
LinearPendulum_plot(T,U,Y,P,F,in,opts,sol)

end


%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 3;

switch num
case 1
    opts = [];
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 600; % number of nodes
case 2
    opts = [];
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
case 3
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50;
    opts.solver.display = 'none';
    opts.solver.tolerance = 1e-15;
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-6;
    opts.dt.meshr.ntmaxinterval = 20;
end

end