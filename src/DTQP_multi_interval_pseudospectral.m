%--------------------------------------------------------------------------
% DTQP_multi_interval_pseudospectral.m
% Implements the multiple-interval pseudospectral method by converting a
% single phase problem into multiphase problem with the appropriate linkage
% constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% TODO: make PS-MI work for multiphase problems (idea: for loop through the
% original phases should work including continuity vs. event constraints as
% needed)
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup,opts] = DTQP_multi_interval_pseudospectral(setup,opts)

% check if this is a usable setup structure
if length(setup) > 1
    error('PS-MI only supports single-phase problems at the moment')
end

% determine the interval boundaries
[t,~,~] = DTQP_MESH_pts(setup(1),opts.dt);

% number of intervals
nphs = length(t);

% update points per interval
opts.dt.nt = opts.dt.nn;

% copy single phase items
setup_single_phase = setup;
UB = setup_single_phase.UB;
LB = setup_single_phase.LB;

% copy single phase setup to all intervals
setup(1:(nphs-1)) = setup_single_phase;

% update options (NOTE: should this handle CGL PS method?)
opts.dt.defects = 'PS';
opts.dt.mesh = 'LGL';

% copy options
opts.dt(1:(nphs-1)) = opts.dt;

%--------------------------------------------------------------------------
% linkage equality constraints
%--------------------------------------------------------------------------
% initialize index
idx = 0;

% determine the number of states (NOTE: this should be done earlier,
% outside this function)
try
   ny = size(setup_single_phase.A,1);
catch
    try
        ny = size(setup_single_phase.B,1);
    catch
        ny = size(setup_single_phase.G,1);
    end
end

% state continuity constraints
e = eye(ny);

% go through each state
for k = 1:ny

    % increment index
    idx = idx + 1;

    % final states
    LY(idx).left.linear.right = 5;
    LY(idx).left.linear.matrix = e(:,k);

    % initial states
    LY(idx).right.linear.right = 4;
    LY(idx).right.linear.matrix = -e(:,k);

    % equality between the final and initial states
    LY(idx).b = 0;

end

% determine the number of controls (NOTE: this should be done earlier,
% outside this function)
try
	nu = size(setup_single_phase.B,2);
catch
    nu = 0;
end

% control continuity constraints (NOTE: expose option?)
e = eye(nu);

% go through each control
for k = 1:nu

    % increment index
    idx = idx + 1;

    % final controls
    LY(idx).left.linear.right = 7;
    LY(idx).left.linear.matrix = e(:,k);

    % initial controls
    LY(idx).right.linear.right = 6;
    LY(idx).right.linear.matrix = -e(:,k);

    % equality between the final and initial controls
    LY(idx).b = 0;

end

% determine the number of parameters (NOTE: this should be done earlier,
% outside this function)
try
	np = size(setup_single_phase.G,2);
catch
    np = 0;
end

% parameter continuity constraints
e = eye(np);

% go through each parameter
for k = 1:np

    % increment index
    idx = idx + 1;

    % final controls
    LY(idx).left.linear.right = 3;
    LY(idx).left.linear.matrix = e(:,k);

    % initial controls
    LY(idx).right.linear.right = 3;
    LY(idx).right.linear.matrix = -e(:,k);

    % equality between the final and initial parameters
    LY(idx).b = 0;

end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% create phases
%--------------------------------------------------------------------------
% combine structures
for phs = 1:nphs-1

    % extract
    UB_right = UB;
    LB_right = LB;

    % remove initial state dependence (4) if not first phase
    if phs ~= 1
        UB_right([UB_right.right] == 4) = [];
        LB_right([LB_right.right] == 4) = [];
    end

    % remove final state dependence (5) if not final phase
    if phs ~= nphs-1
        UB_right([UB_right.right] == 5) = [];
        LB_right([LB_right.right] == 5) = [];
        setup(phs).M = [];
    end

    % assign
    setup(phs).UB = UB_right;
    setup(phs).LB = LB_right;

    % mesh boundary
    setup(phs).t0 = t(phs);
    setup(phs).tf = t(phs+1);

    % linkage constraints
    if phs < nphs-1
        setup(phs).LY = LY;
    end

end
%--------------------------------------------------------------------------

end