%--------------------------------------------------------------------------
% Biegler9p5.m
% Example 9.5 from L. T. Biegler, "Nonlinear Programming." Society for
% Industrial and Applied Mathematics, Jan. 2010.
% doi: 10.1137/1.9780898719383
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

% tunable parameters
tau = 50;
tf = 1; % time horizon

% time horizon
setup.t0 = 0; setup.tf = tf;

% mayer
M(1).left = 5; M(1).right = 5; M(1).matrix = [1,0;0,0];
setup.lq.mayer = M;

% system dynamics
A{1,1} = 0;
A{1,2} = 1;
A{2,1} = tau^2;
A{2,2} = 0;
setup.lq.dynamics.A = A;
setup.lq.dynamics.B = zeros(2,0);
setup.lq.dynamics.Bp = zeros(2,1);
Bv{1,1} = 0;
Bv{2,1} = @(t) - (pi^2+tau^2)*sin(pi*t);
setup.lq.dynamics.Bv = Bv;

% initial states
equality(1).linear(1).right = 4; equality(1).linear(1).matrix = [0;1]; % x0
equality(1).linear(2).right = 3; equality(1).linear(2).matrix = -1; % p
equality(1).b = 0;
equality(2).linear(1).right = 4; equality(2).linear(1).matrix = [1;0]; % x0
equality(2).b = 0;
setup.lq.equality = equality;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Biegler9p5_output(T,U,Y,P,F,in,opts);

% plots
Biegler9p5_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.tolerance = 1e-14;
    opts.solver.display = 'iter';
end

end