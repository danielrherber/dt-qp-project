%--------------------------------------------------------------------------
% Nagurka.m
% Section 6.4 of R. Luus, Iterative Dynamic Programming. CRC Press, 2000,
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


% % auxiliary data
testnum = 2;
switch testnum
    case 1
        n = 6;
    case 2
        n = 20;
    case 3
        n = 50;
end

% time horizon
setup.t0 = 0; setup.tf = 1;


% system dynamics
vec = (1:n);
setup.lq.dynamics.A = [sparse(n-1,1), speye(n-1); ...
    sparse(vec.*(-1).^(vec+1))];
setup.lq.dynamics.B = speye(n);

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(n); % controls
L(2).left = 2; L(2).right = 2; L(2).matrix = eye(n); % states
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5; M(1).matrix = diag([10 zeros(1,n-1)]); % final states
setup.lq.mayer = M;

% simple bounds
UB(1).right = 4; UB(1).matrix = vec';
LB(1).right = 4; LB(1).matrix = vec';
setup.lq.ub = UB; setup.lq.lb = LB;


end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = Nagurka_output(T,U,Y,P,F,in,opts, setup);

% plots
Nagurka_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 200; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
case 2
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 20; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 100;
    opts.solver.display = 'iter';
end

end