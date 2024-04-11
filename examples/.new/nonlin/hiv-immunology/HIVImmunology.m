%--------------------------------------------------------------------------
% HIVImmunology.m
% H. R. Joshi, "Optimal control of an HIV immunology model," Optimal
% Control Applications and Methods, vol. 23, no. 4, pp. 199–213, 2002,
% doi: 10.1002/oca.710
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
tf = 50;
s1 = 2; s2 = 1.5;
m = 0.002; C = 0.007;
k = 2.5e-4; g = 30;
b1 = 14; b2 = 1;
A1 = 2.5e5; A2 = 75;

% time horizon
setup.t0 = 0; setup.tf = tf;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 2;

% system dynamics
str{1} = '[';
str{end+1} = 's1-s2*y2/(b1+y2)- m*y1-k*y1*y2+u1*y1; ';
str{end+1} = 'g*(1-u2)*y2/(b2+y2)-C*y1*y2';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 's1 s2 m C k g b1 b2 A1 A2';
setup.nonlin.data.values = [s1 s2 m C k g b1 b2 A1 A2];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = diag([A1 A2]);
L(2).left = 0; L(2).right = 2; L(2).matrix = -[1 0];
setup.lq.lagrange = L;

% simple bounds
UB(1).right = 4; UB(1).matrix = [400;3];% initial states
LB(1).right = 4; LB(1).matrix = [400;3];
UB(2).right = 2; UB(2).matrix = [1200;5]; % states
LB(2).right = 2; LB(2).matrix = [0;0.05];
UB(3).right = 1; UB(3).matrix = [0.02;0.9]; % controls
LB(3).right = 1; LB(3).matrix = [0;0];
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[400,3];[400,3]];
U0 = [[0.02,0.9];[0,0]];
setup.method.guess.X = [U0,Y0];

% scaling
scaling(1).right = 1; % controls
scaling(1).matrix = [0.02;0.9];
scaling(2).right = 2; % states
scaling(2).matrix = [1200 5];
setup.scaling = scaling;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = HIVImmunology_output(T,U,Y,P,F,in,opts);

% plots
HIVImmunology_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
    case 1
    % default parameters
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
end

end