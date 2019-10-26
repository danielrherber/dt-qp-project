%--------------------------------------------------------------------------
% JadduShimemura.m
% H. Jaddu and E. Shimemura, "Solution of Constrained Linear Quadratic
% Optimal Control Problem Using State Parameterization," Trans. of the
% Society of Instrument and Control Engineers, vol. 34, no. 9, pp.
% 1164-1169, 1998. doi: 10.9746/sicetr1965.34.1164
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = JadduShimemura(varargin)
% input arguments can be provided in the format 'JadduShimemura(p,opts)'

% set local functions
ex_opts = @JadduShimemura_opts; % options function
ex_output = @JadduShimemura_output; % output function
ex_plot = @JadduShimemura_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
p.examplenum = 2; % see below

%% setup
% system dynamics
A = [0,1;0,-1]; 
B = [0;1];

% Lagrange term
L(1).left = 1; % controls
L(1).right = 1; % controls
L(1).matrix = 0.005;
L(2).left = 2; % states
L(2).right = 2; % states
L(2).matrix = eye(2);

% initial states
LB(1).right = 4; % initial states
LB(1).matrix = [0,-1];
UB(1).right = 4; % initial states
UB(1).matrix = [0,-1];

% case number
switch p.examplenum
    case 0 % no additional constraints, not in reference above
        % combine structures
        setup.A = A; setup.B = B; setup.L = L;
        setup.LB = LB; setup.UB = UB; setup.t0 = 0; setup.tf = 1; setup.p = p;
    case 1 % path constraint on state 2
        UB(2).right = 2; % states
        UB(2).matrix = {inf, @(t) 8*(t-0.5).^2 - 0.5};
        
        % combine structures
        setup.A = A; setup.B = B; setup.L = L;
        setup.LB = LB; setup.UB = UB; setup.t0 = 0; setup.tf = 1; setup.p = p;
    case 2 % path constraint on state 1
        UB(2).right = 2; % states
        UB(2).matrix = {@(t) 8*(t-0.5).^2 - 0.5,inf};
        
        % combine structures
        setup.A = A; setup.B = B; setup.L = L;
        setup.LB = LB; setup.UB = UB; setup.t0 = 0; setup.tf = 1; setup.p = p;
    case 3 % intermediate point constraint, multi-phase problem
        LB(2).right = 5; % final states
        LB(2).matrix = [0.5;-inf];
        UB(2).right = 5; % final states
        UB(2).matrix = [0.5;inf];
        
        %--- combine structures, phase 1
        setup(1).A = A; setup(1).B = B; setup(1).L = L;
        setup(1).LB = LB; setup(1).UB = UB; setup(1).t0 = 0; setup(1).tf = 0.5; setup(1).p = p;
        
        %--- combine structures, phase 2 
        setup(2).A = A; setup(2).B = B; setup(2).L = L;
        setup(2).t0 = 0.5; setup(2).tf = 1; setup(2).p = p;
        
        %--- phase 1-2 linkage constraints
        n = length(A);
        clear LY
        q = eye(2);
        for k = 1:n
            LY(k).left.linear.right = 5; % final states
            LY(k).left.linear.matrix = q(:,k);
            LY(k).right.linear.right = 4; % initial states
            LY(k).right.linear.matrix = -q(:,k);
            LY(k).b = 0; % continuous
        end
        setup(1).LY = LY;
end

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = JadduShimemura_opts
% test number
num = 1;

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
end

end