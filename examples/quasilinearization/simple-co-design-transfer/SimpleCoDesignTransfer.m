%--------------------------------------------------------------------------
% SimpleCoDesignTransfer.m
% D. R. Herber, J. T. Allison. 'Nested and simultaneous solution strategies
% for general combined plant and control design problems.' ASME Journal of
% Mechanical Design, 141(1), p. 011402, Jan 2019. doi: 10.1115/1.4040705
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = SimpleCoDesignTransfer(varargin)
% input arguments can be provided in the format 'SimpleCoDesignTransfer(p,opts)'

% set local functions
ex_opts = @SimpleCoDesignTransfer_opts; % options function
ex_output = @SimpleCoDesignTransfer_output; % output function
ex_plot = @SimpleCoDesignTransfer_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
p.x0 = 1; p.v0 = 2; p.tf = 1;

%% setup
% time horizon
p.t0 = 0;

% system dynamics
symb.D = '[y2;-p1*y1 + u1]';
symb.o.ny = 2; % number of states
symb.o.nu = 1; % number of controls
symb.o.np = 1; % number of parameters
symb.o.output = 2; % interp1 compatible

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = 1; % u^2

% simple bounds
UB(1).right = 4; UB(1).matrix = [p.x0;p.v0]; % initial states
LB(1).right = 4; LB(1).matrix = [p.x0;p.v0];
UB(2).right = 5; UB(2).matrix = [0;0]; % final states
LB(2).right = 5; LB(2).matrix = [0;0];

% combine structures
setup.L = L; setup.UB = UB; setup.LB = LB; setup.symb = symb;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
% disp("paused"); pause
% ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = SimpleCoDesignTransfer_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts.general.displevel = 1;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 40; % number of nodes
case 2
    opts.general.displevel = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 40; % number of nodes
case 3
    opts.general.displevel = 1;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 40; % number of nodes
    opts.qlin.sqpflag = false;
    opts.qlin.trustregionflag = true;
    opts.qlin.delta = 0.25;
end

end