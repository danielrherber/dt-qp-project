%--------------------------------------------------------------------------
% TetherAssistedReentry.m
% Example 1 in P. Williams, "Comparison of Differentiation and Integration
% Based Direct Transcription Methods." AAS/AIAA Space Flight Mechanics
% Meeting, AAS 05-128, Jan. 2005
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
setup.t0 = 0; setup.tf = 5;

% counts for controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 4;
setup.counts.np = 0;
setup.counts.nv = 0;

% lagrange
setup.nonlin.lagrange = '(y3*( (y2+1)^2 + 3*cos(y1)^2 - 1 ) - u1)^2/2';

% dynamics
strD{1} = '[';
strD{end+1} = 'y2;';
strD{end+1} = '-2*(y2+1)*y4/y3-3*sin(y1)*cos(y1);';
strD{end+1} = 'y4;';
strD{end+1} = 'y3*((y2+1)^2 + 3*cos(y1)^2 - 1) - u1';
strD{end+1} = ']';
setup.nonlin.dynamics = horzcat(strD{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = [deg2rad(30),0,0.05,0.05]; % initial states
LB(1).right = 4; LB(1).matrix = [deg2rad(30),0,0.05,0.05];
UB(2).right = 5; UB(2).matrix = [0,sqrt(6)/2,1,0]; % final states
LB(2).right = 5; LB(2).matrix = [0,sqrt(6)/2,1,0];
UB(3).right = 1; UB(3).matrix = 12; % controls
LB(3).right = 1; LB(3).matrix = 0.01;
LB(4).right = 2; LB(4).matrix = [-inf,-inf,-inf,0]; % states
setup.lq.ub = UB; setup.lq.lb = LB;

% guess
Y0 = [[deg2rad(30),0,0.05,0.05];[0,sqrt(6),1,0]];
U0 = [0.01;0.01];
P0 = [];
setup.method.guess.X = [U0,Y0,P0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = TetherAssistedReentry_output(T,U,Y,P,F,in,opts);

% plots
TetherAssistedReentry_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
    opts.solver.tolerance = 1e-12;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50;
    opts.solver.tolerance = 1e-12;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
end

end
