%--------------------------------------------------------------------------
% MinimumEnergyTransfer.m
% pp. 21-23 of J.  Klamka, Controllability and Minimum Energy Control,
% Springer, 2019, doi: 10.1007/978-3-319-92540-0
% pp. 163-164 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
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
ny = 18; % number of states
nu = 6; % number of controls
rng(83233683,'twister') % random number seed
Adensity = rand;
Aeig = -2 + (2 - -2).*rand(ny,1);
auxdata.A = sprandsym(ny,Adensity,Aeig);
auxdata.B = -10 + (10 - -10).*rand(ny,nu);
auxdata.y0 = 10*rand(ny,1); % initial states
auxdata.yf = zeros(ny,1); % final states

setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 2;

% system dynamics
setup.lq.dynamics.A = auxdata.A;
setup.lq.dynamics.B = auxdata.B;

% check controllability
try
	Co = ctrb(auxdata.A,auxdata.B);
    if rank(Co) ~= ny
        warning('system is not controllable')
    end
catch
    warning('unable to check controllability')
end

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = eye(nu); % control variables
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4;  LB(1).matrix = auxdata.y0; % initial states
UB(1).right = 4;  UB(1).matrix = auxdata.y0; % initial states
LB(2).right = 5;  LB(2).matrix = auxdata.yf; % final states
UB(2).right = 5;  UB(2).matrix = auxdata.yf; % final states
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = MinimumEnergyTransfer_output(T,U,Y,P,F,in,opts);

% plots
MinimumEnergyTransfer_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
case 2
    opts = [];
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50; % number of nodes
case 3
    opts = [];
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 4; % number of nodes
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-8;
end

end