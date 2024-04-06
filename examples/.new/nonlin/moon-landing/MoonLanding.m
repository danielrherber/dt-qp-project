%--------------------------------------------------------------------------
% MoonLanding.m
% Example 2 in P. Williams, "Comparison of Differentiation and Integration
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
setup.t0 = 0;
setup.tf = 1;

% counts for controls, states, and parameters
setup.counts.nu = 1;
setup.counts.nx = 3;
setup.counts.np = 1;

% mayer
M(1).left = 0; % singleton
M(1).right = 5; M(1).matrix = -[0,0,1];
setup.lq.mayer = M;

% dynamics
strD{1} = '[';
strD{end+1} = 'p1*(y2);';
strD{end+1} = 'p1*(-1+u1/y3);';
strD{end+1} = 'p1*(-u1/2.3)';
strD{end+1} = ']';
setup.nonlin.dynamics = horzcat(strD{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,-0.783,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,-0.783,1];
UB(2).right = 5; UB(2).matrix = [0,0,1]; % final states
LB(2).right = 5; LB(2).matrix = [0,0,0];
UB(3).right = 3; UB(3).matrix = 10; % parameters
LB(3).right = 3; LB(3).matrix = 0;
UB(3).right = 1; UB(3).matrix = 1.1; % controls
LB(3).right = 1; LB(3).matrix = 0;
setup.lq.ub = UB; setup.lq.lb = LB;

% initial guess
Y0 = [[1,-0.783,1];[0,0,1]];
U0 = [1;0];
P0 = [0.5;0.5];
setup.method.guess.X = [U0,Y0,P0];

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = MoonLanding_output(T,U,Y,P,F,in,opts);

% plots
MoonLanding_plot(T,U,Y,P,F,in,opts,sol)

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
    opts.dt.nt = 1000;
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