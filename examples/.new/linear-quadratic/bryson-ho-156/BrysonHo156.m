%--------------------------------------------------------------------------
% BrysonHo156.m
% p. 156 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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
[O,sol] = postProcessing(T,U,Y,P,F,in,opts);


%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% auxiliary data
auxdata.x0 = 10; auxdata.v0 = 1; 
auxdata.c = 1; auxdata.omega = 1;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 10;

% system dynamics
setup.lq.dynamics.A = [0,1;-auxdata.omega^2,0];
setup.lq.dynamics.B = [0;1];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = [auxdata.c/2,0;0,0];
setup.lq.mayer = M;

% simple bounds
LB(1).right = 4; LB(1).matrix = [auxdata.x0;auxdata.v0]; % initial states
UB(1).right = 4; UB(1).matrix = [auxdata.x0;auxdata.v0];
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = BrysonHo156_output(T,U,Y,P,F,in,opts);

% plots
BrysonHo156_plot(T,U,Y,P,F,in,opts,sol)

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
    opts = [];
case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 5;
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-8;
end

end