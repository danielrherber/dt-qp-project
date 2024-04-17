%--------------------------------------------------------------------------
% DTQP3.m
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
auxdata.x0 = 1;  auxdata.b = 1;  auxdata.r = 1;  auxdata.m = 1;
auxdata.a1 = 2;  auxdata.a2 = 1; auxdata.w1 = 3; auxdata.w2 = 8;
auxdata.ParameterFlag = 1; % parameter version?
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 10;

% system dynamics
if auxdata.ParameterFlag
	setup.lq.dynamics.Bp = auxdata.b;
else
	setup.lq.dynamics.B = auxdata.b;
end
setup.lq.dynamics.Bv{1,1} = @(t) auxdata.a1*sin(auxdata.w1*t) + auxdata.a2*sin(auxdata.w2*t);


% Lagrange term
if auxdata.ParameterFlag
    L(1).left = 3; L(1).right = 3;  L(1).matrix = auxdata.r; % parameters
else
    L(1).left = 1; L(1).right = 1; L(1).matrix = auxdata.r; % control variables
end
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; M(1).right = 5;  M(1).matrix = auxdata.m; % final states
setup.lq.mayer = M;

% simple bounds
LB(1).right = 4; LB(1).matrix = auxdata.x0;  % initial states
UB(1).right = 4; UB(1).matrix = auxdata.x0; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end


%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = DTQP3_output(T,U,Y,P,F,in,opts);

% plots
DTQP3_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 100;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
end

end