%--------------------------------------------------------------------------
% JadduShimemura.m
% H. Jaddu and E. Shimemura, "Solution of Constrained Linear Quadratic
% Optimal Control Problem Using State Parameterization," Trans. of the
% Society of Instrument and Control Engineers, vol. 34, no. 9, pp.
% 1164-1169, 1998. doi: 10.9746/sicetr1965.34.1164
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
auxdata.examplenum = 3; % see below
setup.auxdata = auxdata;

% system dynamics
setup.lq.dynamics.A = [0,1;0,-1];
setup.lq.dynamics.B = [0;1];


% Lagrange term
L(1).left = 1; % controls
L(1).right = 1; % controls
L(1).matrix = 0.005;
L(2).left = 2; % states
L(2).right = 2; % states
L(2).matrix = eye(2);
setup.lq.lagrange = L;

% simple bounds
LB(1).right = 4;  LB(1).matrix = [0,-1]; % initial states
UB(1).right = 4;  UB(1).matrix = [0,-1]; % initial states


% case number
switch auxdata.examplenum
    case 0 % no additional constraints, not in reference above
        
        % time horizon
        setup.t0 = 0; setup.tf = 1; 

        setup.lq.ub = UB; setup.lq.lb = LB;

    case 1 % path constraint on state 2
        
        % time horizon
        setup.t0 = 0; setup.tf = 1; 

        UB(2).right = 2; % states
        UB(2).matrix = {inf, @(t) 8*(t-0.5).^2 - 0.5};
        setup.lq.ub = UB; setup.lq.lb = LB;
        
    case 2 % path constraint on state 1

        % time horizon
        setup.t0 = 0; setup.tf = 1;

        UB(2).right = 2; % states
        UB(2).matrix = {@(t) 8*(t-0.5).^2 - 0.5,inf};
        setup.lq.ub = UB; setup.lq.lb = LB;

    case 3 % intermediate point constraint, multi-phase problem

        % initialize struct
        setup_ = DTQP_setup_initialize;
        setup(2) = setup_;

        % time horizon
        setup(1).t0 = 0; setup(1).tf = 0.5;

        LB(2).right = 5; % final states
        LB(2).matrix = [0.5;-inf];
        UB(2).right = 5; % final states
        UB(2).matrix = [0.5;inf];
        setup(1).lq.ub = UB; setup(1).lq.lb = LB;

        %--- combine structures, phase 2
        setup(2).lq.dynamics.A = setup(1).lq.dynamics.A; 
        setup(2).lq.dynamics.B = setup(1).lq.dynamics.B;
        setup(2).lq.lagrange.L = L;
        setup(2).t0 = 0.5; setup(2).tf = 1; 
        setup(2).auxdata = auxdata;

        %--- phase 1-2 linkage constraints
        % n = length(setup(1).lq.dynamics.A);
        % clear LY
        % q = eye(2);
        % for k = 1:n
        %     LY(k).left.linear.right = 5; % final states
        %     LY(k).left.linear.matrix = q(:,k);
        %     LY(k).right.linear.right = 4; % initial states
        %     LY(k).right.linear.matrix = -q(:,k);
        %     LY(k).b = 0; % continuous
        % end
        % setup(1).lq.linkage_equality = LY;
end

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% outputs
[O,sol] = JadduShimemura_output(T,U,Y,P,F,in,opts);

% plots
JadduShimemura_plot(T,U,Y,P,F,in,opts,sol)

end



%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 2;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 50;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000;
case 4
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
    opts.solver.tolerance = 1e-15;
    opts.solver.display = 'none';
    opts.dt.meshr.method = 'SS-BETTS';
    opts.dt.meshr.tolerance = 1e-4;
end

end