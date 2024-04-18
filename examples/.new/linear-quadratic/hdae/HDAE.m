%--------------------------------------------------------------------------
% HDAE.m
% Heat Diffusion Process with Inequality
%--------------------------------------------------------------------------
% pp. 192-195 of J. T. Betts, Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming. SIAM, 2010, isbn: 9780898716887.
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


%% tunable parameters
n = 20;    a = 0.5;    b = 0.2;    c = 1;
q1 = 1e-3; q2 = 1e-3;  delta = pi/n;

% auxiliary data
auxdata.n = n;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 5;

% system dynamics
A = gallery('tridiag',n-1,1,-2,1);
setup.lq.dynamics.A = A/delta^2;
B = zeros(n-1,2);
B(1,1) = 1;
B(n-1,2) = 1;
setup.lq.dynamics.B = B/delta^2;

% Lagrange term
idx = 0;

% quadratic control
idx = idx + 1;
L(idx).left = 1; % controls;
L(idx).right = 1; % controls;
L(idx).matrix = diag([1/2*delta+q1,1/2*delta+q2]);

% quadratic states
idx = idx + 1;
L(idx).left = 2; % states
L(idx).right = 2; % states
L(idx).matrix = delta*eye(n-1);
setup.lq.lagrange = L;


% simple bounds
idx = 0;

% initial states
idx = idx + 1;
LB(idx).right = 4; % initial states
LB(idx).matrix = zeros(n-1,1);

% state inequality constraints
idx = idx + 1;
LB(idx).right = 2; % states
g = cell(1,n-1);
for k = 1:n-1
   g{1,k} = @(t)  c*(sin(k*pi/n)*sin(pi*t/5) - a) - b;
end
LB(idx).matrix = g;

% control inequality constraints
idx = idx + 1;
LB(idx).right = 1; % controls
LB(idx).matrix{1,1} = @(t)  c*(sin(0*pi/n)*sin(pi*t/5) - a) - b;
LB(idx).matrix{2,1} = @(t)  c*(sin(n*pi/n)*sin(pi*t/5) - a) - b;

% initial controls
idx = idx + 1;
LB(idx).right = 1; % controls
LB(idx).matrix{1,1} = @(t) [0;-inf(length(t)-1,1)];
LB(idx).matrix{2,1} = @(t) [0;-inf(length(t)-1,1)];

%--- simple upper bounds
idx = 0;

% initial states
idx = idx + 1;
UB(idx).right = 4; % initial states
UB(idx).matrix = zeros(n-1,1);

% initial controls
idx = idx + 1;
UB(idx).right = 1; % controls
UB(idx).matrix{1,1} = @(t) [0;inf(length(t)-1,1)];
UB(idx).matrix{2,1} = @(t) [0;inf(length(t)-1,1)];
setup.lq.ub = UB; 
setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = HDAE_output(T,U,Y,P,F,in,opts);

% plots
HDAE_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 100;
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
case 3
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-2;
end

end