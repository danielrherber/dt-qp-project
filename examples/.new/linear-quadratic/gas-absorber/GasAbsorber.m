%--------------------------------------------------------------------------
% GasAbsorber.m
% Section 7.4.4 of R. Luus, Iterative Dynamic Programming. CRC Press, 2000,
% isbn: 1584881488
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
testnum = 3;
switch testnum
    case 1
    n = 10;
    case 2
    n = 25;
    case 3
    n = 50;
    case 4
    n = 100;
    case 5
    n = 200;
end

% time horizon
setup.t0 = 0; setup.tf = 15;

% system dynamics
setup.lq.dynamics.A = spdiags([0.538998*ones(n,1) -1.173113*ones(n,1) 0.634115*ones(n,1)],-1:1,n,n);
setup.lq.dynamics.B = sparse(n,2);
setup.lq.dynamics.B(1,1)   = 0.538998;
setup.lq.dynamics.B(end,2) = 0.634115;

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(2); % controls
L(2).left = 2; L(2).right = 2; L(2).matrix = eye(n); % states
setup.lq.lagrange = L;

% simple bounds
y0 = -0.0307-(0.1273-0.0307)/(n-1)*((1:n)'-1);
UB(1).right = 4; UB(1).matrix = y0;  % initial states
LB(1).right = 4; LB(1).matrix = y0;  % initial states
LB(2).right = 1; LB(2).matrix = [0 0]; % control bounds
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = GasAbsorber_output(T,U,Y,P,F,in,setup,opts);

% plots
GasAbsorber_plot(T,U,Y,P,F,in,opts,sol)

end



%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 1;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
end

end