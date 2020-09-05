%--------------------------------------------------------------------------
% TransferMaxRadius.m
% p. 66-69 of A. E. Bryson Jr. and Y.-C. Ho, Applied Optimal Control.
% Taylor & Francis, 1975, isbn: 9780891162285
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = TransferMaxRadius(varargin)
% input arguments can be provided in the format 'TransferMaxRadius(p,opts)'

% set local functions
ex_opts = @TransferMaxRadius_opts; % options function
ex_output = @TransferMaxRadius_output; % output function
ex_plot = @TransferMaxRadius_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
xT  = 0.1405;
m0 = 1;
dm = 0.0749;
xmu = 1;
tf = 3.32;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 2; n.ny = 4;

% system dynamics
str{1} = '[';
str{end+1} = 'y3; ';
str{end+1} = 'y4/y1; ';
str{end+1} = 'y4^2/y1 - xmu/y1^2 + (xT/(m0 - dm*t))*u1; ';
str{end+1} = '-y3*y4/y1 + (xT/(m0 - dm*t))*u2';
str{end+1} = ']';
symb.D = horzcat(str{:});
symb.paramstr = 'xmu xT m0 dm';
symb.param = [xmu xT m0 dm];

% simple bounds
UB(1).right = 4; UB(1).matrix = [1,0,0,1]; % initial states
LB(1).right = 4; LB(1).matrix = [1,0,0,1];
UB(2).right = 5; UB(2).matrix = [inf,inf,0,inf]; % final states
LB(2).right = 5; LB(2).matrix = [-inf,-inf,0,-inf];

% equality constraints
symb.ceq.func = '[u1^2 + u2^2 - 1; sqrt(xmu/yf1) - yf4]';
symb.ceq.pathboundary = [1 0];

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-1,0,0,0];

% guess
Y0 = [[1,0,0,1];[1.5,pi,0,0.5]];
U0 = [[0,1];[-1,0]];
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.UB = UB; setup.LB = LB; setup.M = M;
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
function opts = TransferMaxRadius_opts
% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 40; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.method.form = 'nonlinearprogram';
end

end