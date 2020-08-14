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

casenum = 2;

switch casenum
    case 1
        % system dynamics
        symb.D = '[y2;-y1+y2-y1^2*y2+u1]';
        n.ny = 2; % number of states
        n.nu = 1; % number of inputs

        % Lagrange term
        symb.Ob = 'y1^2 + y2^2 + u1^2';
        % L(1).left = 1; L(1).right = 1; L(1).matrix = 1;
        % L(2).left = 2; L(2).right = 2; L(2).matrix = eye(2);

        % simple bounds
        UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
        LB(1).right = 4; LB(1).matrix = [1;0];
        UB(2).right = 1; UB(2).matrix = 1; % controls
        LB(2).right = 1; LB(2).matrix = -0.3;

        % guess
        Y0 = [[1,0];[1,0]];
        U0 = [[0];[0]];
        p.guess = [U0,Y0];

        % combine structures
        setup.symb = symb; setup.UB = UB; setup.LB = LB; % setup.L = L;
        setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

    case 2 % co-design example
        dmin = [0.1 0.1];
        dmax = [5 5];

        % system dynamics
        symb.D = '[y2;-p1*y1+y2-p2*y1^2*y2+u1]';
        n.ny = 2; % number of states
        n.nu = 1; % number of inputs
        n.np = 2;

        % Lagrange term
        if ~isfield(opts.method,'olqflag') || opts.method.olqflag
            L(1).left = 1; L(1).right = 1; L(1).matrix = 1;
            L(2).left = 2; L(2).right = 2; L(2).matrix = eye(2);
            setup.L = L;
        else
            symb.Ob = 'y1^2 + y2^2 + u1^2';
        end

        % simple bounds
        UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
        LB(1).right = 4; LB(1).matrix = [1;0]; % initial states
        UB(2).right = 1; UB(2).matrix = 1; % controls
        LB(2).right = 1; LB(2).matrix = -0.5; % controls
        UB(3).right = 3; UB(3).matrix = dmax; % parameters
        LB(3).right = 3; LB(3).matrix = dmin; % parameters
        LB(4).right = 2; LB(4).matrix = [-inf;-0.4]; % states

        % guess
        Y0 = [[1,0];[0,0]];
        U0 = [[0];[0]];
        P0 = [[1,1];[1,1]];
        p.guess = [U0,Y0,P0];

        % combine structures
        setup.symb = symb; setup.UB = UB; setup.LB = LB;
        setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

    case 3 % completely symbolic formulation
        % system dynamics
        symb.D = '[y2; -y1 + y2 - y1^2*y2 + u1]';
        n.ny = 2; % number of states
        n.nu = 1; % number of inputs

        % Lagrange term
        symb.Ob = 'y1^2 + y2^2 + u1^2';

        % equality constraints
        symb.ceq = '[yi1 - 1; yi2 - 0]';

        % inequality constraints
        symb.cin = '[u1 - 1; -u1 - 0.3]';

        % guess
        Y0 = [[1,0];[1,0]];
        U0 = [[0];[0]];
        p.guess = [U0,Y0];

        % combine structures
        setup.symb = symb; % setup.UB = UB; setup.LB = LB; % setup.L = L;
        setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;
end

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
% disp("paused"); pause % for quasilinearization plots
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = Vanderpol_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.method.form = 'nonlinearprogram';
    opts.solver.function = 'ipfmincon';
    opts.solver.maxiters = 4000;
    opts.solver.display = 'iter';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 40; % number of nodes
    opts.method.form = 'nonlinearprogram';
    opts.solver.function = 'ipfmincon';
    opts.solver.display = 'iter';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
case 3
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.method.maxiters = 5000;
    opts.dt.nt = 40; % number of nodes
    opts.method.form = 'qlin';
    opts.solver.display = 'none';
    opts.method.trustregionflag = false;
    opts.method.improveguess = false;
end

end