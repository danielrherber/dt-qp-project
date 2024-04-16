%--------------------------------------------------------------------------
% DTQP1.m
% Complex linear-quadratic dynamic optimization problem
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
auxdata.g = @(t) sin(2*pi*t) + 0.5;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 1;

% system dynamics
setup.lq.dynamics.A = [-1,2,0,0;3,-4,0,0;1,2,-1,0;1,0,0,0];
setup.lq.dynamics.B = [1,0;-1,0;0,1/20;0,0];
setup.lq.dynamics.Bp = zeros(4,1);

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(2)/10; % u1^2 + u2^2
L(2).left = 1; L(2).right = 2; L(2).matrix = [1,1,0,0;0,0,0,0]; % u1*y1 + u1*y2
L(3).left = 2; L(3).right = 2; L(3).matrix = zeros(4); L(3).matrix(2,2) = 5; % 5*y2^2
L(4).left = 0; L(4).right = 2; L(4).matrix = {0,@(t) -5*2*auxdata.g(t),0,0}; % -5*2*g*y2
L(5).left = 0; L(5).right = 0; L(5).matrix{1} = @(t) 5*(auxdata.g(t)).^2; % 5*g^2
setup.lq.lagrange = L;

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = 1; % p1
setup.lq.mayer = M;

% equality constraint: y2(t0)-y2(tf) = 0, equality constraint
Y(1).linear(1).right = 4; Y(1).linear(1).matrix = [0;1;0;0]; % y2(t0)
Y(1).linear(2).right = 5; Y(1).linear(2).matrix = [0;-1;0;0]; % -y2(tf)
Y(1).b = 0;
setup.lq.equality = Y;

% inequality constraint: -y1 + u2/12 < 0
Z(1).linear(1).right = 2; Z(1).linear(1).matrix = [-1;0;0;0]; % -y1
Z(1).linear(2).right = 1; Z(1).linear(2).matrix = [0;1/12]; % u2/12
Z(1).b = 0;

% inequality constraint: y3 < p
Z(2).linear(1).right = 2; Z(2).linear(1).matrix = [0;0;1;0]; % y3
Z(2).linear(2).right = 3; Z(2).linear(2).matrix = -1; % -p1
Z(2).b = 0;

setup.lq.inequality = Z;

% simple bounds
UB(1).right = 4; UB(1).matrix = [2;inf;0.5;0];    % initial states
LB(1).right = 4; LB(1).matrix = [2;-inf;0.5;0];   % initial states
UB(2).right = 5; UB(2).matrix = [inf;inf;inf;0];  % final states
LB(2).right = 5; LB(2).matrix = -[inf;inf;inf;0]; % final states
UB(3).right = 1; UB(3).matrix = [inf;10];  % abs(u2) < 10, simple bounds
LB(3).right = 1; LB(3).matrix = -[inf;10]; % abs(u2) < 10, simple bounds
UB(4).right = 2; UB(4).matrix= {inf;@(t) auxdata.g(t);inf;inf}; % y2 < g(t)
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = DTQP1_output(T,U,Y,P,F,in,opts);

% plots
DTQP1_plot(T,U,Y,P,F,in,opts,sol)

end



%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 2;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 5000; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 20; % number of nodes
    opts.method.reordervariables = 0;
    opts.solver.function = 'built-in';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiters = 100;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-5;
end

end