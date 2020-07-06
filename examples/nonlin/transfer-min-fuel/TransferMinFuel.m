%--------------------------------------------------------------------------
% TransferMinFuel.m
% I. M. Ross, Q. Gong, and P. Sekhavat, "Low-Thrust, High-Accuracy
% Trajectory Optimization," Journal of Guidance, Control, and Dynamics,
% vol. 30, no. 4, pp. 921–933, Jul. 2007, doi: 10.2514/1.23181
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = TransferMinFuel(varargin)
% input arguments can be provided in the format 'TransferMinFuel(p,opts)'

% set local functions
ex_opts = @TransferMinFuel_opts; % options function
ex_output = @TransferMinFuel_output; % output function
ex_plot = @TransferMinFuel_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
testnum = 1; % see below
umax = 0.01;

%% setup
% time horizon
p.t0 = 0; p.tf = 57;

% number of controls, states, and parameters
n.nu = 2; n.ny = 4;

% system dynamics
str{1} = '[';
str{end+1} = 'y3; ';
str{end+1} = 'y4/y1; ';
str{end+1} = 'y4^2/y1 - 1/y1^2 + u1; ';
str{end+1} = '-y3*y4/y1 + u2';
str{end+1} = ']';
symb.D = horzcat(str{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,0,0,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,0,0,1];
UB(2).right = 5; UB(2).matrix = [4,inf,0,0.5]; % final states
LB(2).right = 5; LB(2).matrix = [4,-inf,0,0.5];
UB(3).right = 2; UB(3).matrix = [4,inf,0.5,1]; % states
LB(3).right = 2; LB(3).matrix = [1,-inf,0,0];

switch testnum
    %----------------------------------------------------------------------
    case 1 % p = 1, q = inf

    % Lagrange term
    % symb.Ob = 'abs(u1) + abs(u2)';
    symb.Ob = 'sqrt(u1^2) + sqrt(u2^2)'; % works with complex numbers

    % linear inequality constraint
    UB(4).right = 1; UB(4).matrix = [umax, umax]; % controls
    LB(4).right = 1; LB(4).matrix = [-umax, -umax];
    %----------------------------------------------------------------------
    case 2 % p = 2, q = 2

    % Lagrange term
    symb.Ob = '(u1^2 + u2^2)^(1/2)';
    % symb.Ob = 'u1^2 + u2^2';
    % L(1).left = 1; L(1).right = 1; L(1).matrix = eye(2); setup.L = L; % controls

    % nonlinear inequality constraint
    symb.cin.func = 'u1^2 + u2^2 - umax^2';
    symb.cin.pathboundary = 1;
    symb.paramstr = 'umax';
    symb.param = umax;

    % linear inequality constraint
    UB(4).right = 1; UB(4).matrix = [umax, umax]; % controls
    LB(4).right = 1; LB(4).matrix = [-umax, -umax];
end

% guess
Y0 = [[1,0,0,1];[4,0,0,0.5]];
U0 = [[0.01,0.01];[0.01,0.01]];
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.UB = UB; setup.LB = LB;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

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
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = TransferMinFuel_opts
% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.tolerance = 1e-10;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end