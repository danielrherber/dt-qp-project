%--------------------------------------------------------------------------
% Hager1.m
% W. W. Hager, Runge-Kutta Methods in Optimal Control and the Transformed
% Adjoint System," Numerische Mathematik, vol. 87, no. 2, pp. 247-282,
% Dec. 2000. doi: 10.1007/s002110000178
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

% time horizon
setup.t0 = 0; setup.tf = 1;

% system dynamics
setup.lq.dynamics.A = 0.5;
setup.lq.dynamics.B = 1;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1; 
L(2).left = 2; % state variables
L(2).right = 1; % control variables
L(2).matrix = 1;
L(3).left = 2; % state variables
L(3).right = 2; % state variables
L(3).matrix = 5/4;
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4;  LB(1).matrix = 1; % initial states
UB(1).right = 4;  UB(1).matrix = 1; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;


end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = Hager1_output(T,U,Y,P,F,in,opts);

% plots
Hager1_plot(T,U,Y,P,F,in,opts,sol)

end


%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
end

end