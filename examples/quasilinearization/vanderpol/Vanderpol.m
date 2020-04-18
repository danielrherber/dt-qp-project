%--------------------------------------------------------------------------
% Vanderpol.m
% E. B. Canto, et al., "Restricted second order information for the
% solution of optimal control problems using control vector
% parameterization", *Journal of Process Control*, vol. 12, no. 2002,
% pp. 243-255, 2002, doi: 10.1016/S0959-1524(01)00008-7
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Vanderpol(varargin)
% input arguments can be provided in the format 'Vanderpol(p,opts)'

% set local functions
ex_opts = @Vanderpol_opts; % options function
ex_output = @Vanderpol_output; % output function
ex_plot = @Vanderpol_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon
p.t0 = 0; p.tf = 5;

casenum = 1;

switch casenum
    case 1
        % system dynamics
        symb.D = '[y2;-y1+y2-y1^2*y2+u1]';
        symb.o.ny = 2; % number of states
        symb.o.nu = 1; % number of inputs
        symb.o.output = 2; % interp1 compatible

        % Lagrange term
        symb.Ob = 'y1^2 + y2^2 + u1^2';

        % simple bounds
        UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
        LB(1).right = 4; LB(1).matrix = [1;0];
        UB(2).right = 1; UB(2).matrix = 1; % controls
        LB(2).right = 1; LB(2).matrix = -0.3;

        % combine structures
        setup.symb = symb; setup.UB = UB; setup.LB = LB;
        setup.t0 = p.t0; setup.tf = p.tf; setup.p = p;
    case 2 % co-design example
        dmin = [0.1 0.1];
        dmax = [5 5];
        % dmin = [0.153 1.54]*1;
        % dmax = [0.153 1.54]*1;

        % system dynamics
        symb.D = '[y2;-p1*y1+y2-p2*y1^2*y2+u1]';
        symb.o.ny = 2; % number of states
        symb.o.nu = 1; % number of inputs
        symb.o.np = 2; % number of parameters
        symb.o.output = 2; % interp1 compatible

        % Lagrange term
        symb.Ob = 'y1^2 + y2^2 + u1^2';

        % simple bounds
        UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
        LB(1).right = 4; LB(1).matrix = [1;0]; % initial states
        UB(2).right = 1; UB(2).matrix = 1; % controls
        LB(2).right = 1; LB(2).matrix = -0.5; % controls
        UB(3).right = 3; UB(3).matrix = dmax; % parameters
        LB(3).right = 3; LB(3).matrix = dmin; % parameters
        LB(4).right = 2; LB(4).matrix = [-0.4;-inf]; % states

        % combine structures
        setup.symb = symb; setup.UB = UB; setup.LB = LB;
        setup.t0 = p.t0; setup.tf = p.tf; setup.p = p;

end

%% solve
t1 = tic;
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);
toc(t1);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
% disp("paused"); pause % for quasilinearization plots
% ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = Vanderpol_opts
% test number
num = 3;

switch num
case 1
    % default parameters
    opts.general.displevel = 1;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 400; % number of nodes
case 2
    opts.general.displevel = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 40; % number of nodes
    opts.qlin.sqpflag = false;
    opts.qlin.imax = 20;
case 3
    % default parameters
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.qlin.sqpflag = false;
    opts.qlin.imax = 20;
end

end