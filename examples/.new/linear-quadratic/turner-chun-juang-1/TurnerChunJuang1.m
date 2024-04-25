%--------------------------------------------------------------------------
% TurnerChunJuang1.m
% J. D. Turner, H. M. Chun, J. N. Juang, "Closed-Form Solutions for a Class
% of Optimal Quadratic Tracking Problems", Journal of Optimization Theory
% and Applications, vol. 47, no. 4, pp. 465-481, Dec. 1985.
% doi: 10.1007/BF00942192
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
auxdata.y0 = -5; % initial state
auxdata.em = 1;
auxdata.um = 1;
auxdata.a = 1;
auxdata.S = 5;
auxdata.eta = 10;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = 10;

% system dynamics
setup.lq.dynamics.A = -1/auxdata.a;
setup.lq.dynamics.B = 1;

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = auxdata.um^-2;
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = auxdata.em^-2;
L(3).left = 0; % singleton
L(3).right = 2; % state variables
L(3).matrix = {@(t) -2*auxdata.em^-2*auxdata.eta};
L(4).left = 0; % singleton
L(4).right = 0; % singleton
L(4).matrix = {@(t) auxdata.em^-2*auxdata.eta.^2};
setup.lq.lagrange = L;

% Mayer term
M(1).left = 5; % final states
M(1).right = 5; % final states
M(1).matrix = auxdata.S;
M(2).left = 0; % singleton
M(2).right = 5; % final states
M(2).matrix = -2*auxdata.S*auxdata.eta;
M(3).left = 0; % singleton
M(3).right = 0; % singleton
M(3).matrix = auxdata.S*auxdata.eta^2;
setup.lq.mayer = M;

% simple bounds
LB(1).right = 4;  LB(1).matrix = auxdata.y0; % initial states
UB(1).right = 4;  UB(1).matrix = auxdata.y0; % initial states
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = TurnerChunJuang1_output(T,U,Y,P,F,in,opts);

% plots
TurnerChunJuang1_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 20;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'C';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
end

end