%--------------------------------------------------------------------------
% SimpleSuspensionSimultaneous.m
% D. R. Herber and A. K. Sundarrajan, "On the uses of linear-quadratic
% methods in solving nonlinear dynamic optimization problems with direct
% transcription", in ASME International Mechanical Engineering Congress &
% Exposition, 2020.
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = SimpleSuspensionSimultaneous(varargin)

% set local functions
ex_opts = @SimpleSuspensionSimultaneous_opts;
ex_output = @SimpleSuspension_output;
ex_plot = @SimpleSuspension_plot;

% set p and opts
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

% problem parameters
p = SimpleSuspensionProblemParameters;

% tunable parameters
p.t0 = 0; p.tf = 3;

%% setup
% number of controls, states, and parameters
n.nu = 1; n.ny = 4; n.np = 2;

% system dynamics
str{1} =  '[';
str{end+1} = 'y2 - z0d;';
str{end+1} = '(- kt*y1 - (p1 + ct)*y2 + p2*y3 + p1*y4 + ct*z0d - u1)/mus;'; % missing y2 term
str{end+1} = 'y4 - y2;';
str{end+1} = '( p1*y2 - p2*y3 - p1*y4 + u1)/ms';
str{end+1} = ']';
element.dynamics = horzcat(str{:});

% objective function
str0{1} = 'w1*y1^2 + ';
str0{2} = 'w2/ms^2*( p1*y2 - p2*y3 - p1*y4 + u1 )^2 + ';
str0{3} = 'w3*u1^2';
element.lagrange = horzcat(str0{:});

% symbolic parameters
element.parameter_list = 'ct kt mus ms w1 w2 w3 z0d';
element.parameter_values = {p.bt p.kt p.mu p.ms p.w1 p.w2 p.w3 p.z0dot};

% initial state values
LB(1).right = 4;
LB(1).matrix = zeros(4,1);
UB(1).right = 4;
UB(1).matrix = zeros(4,1);

% simple parameter bounds
LB(2).right = 3;
LB(2).matrix = [p.bmin;p.kmin];
UB(2).right = 3;
UB(2).matrix = [p.bmax;p.kmax];

% rattlespace constraints
LB(3).right = 2; % states
LB(3).matrix = [-inf,-inf,-p.rmax,-inf];
UB(3).right = 2; % states
UB(3).matrix = [inf,inf,p.rmax,inf];

% guess
Y0 = [[0,0,0,0];[0,0,0,0]];
U0 = [[0];[0]];
P0 = [[LB(2).matrix'];[LB(2).matrix']];
setup.guess.X = [U0,Y0,P0];

% scaling (with knowledge of the solution)
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = 400;
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [0.004 0.4 p.rmax 0.04];
setup.scaling(3).right = 3; % parameters
setup.scaling(3).matrix = [1e2;1e4];

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

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
function opts = SimpleSuspensionSimultaneous_opts

opts.general.displevel = 2;
opts.general.plotflag = 1;
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.dt.nt = 1000;
opts.method.form = 'nonlinearprogram';
opts.method.olqflag = true;
opts.solver.tolerance = 1e-14;
opts.solver.maxiters = 10000;
opts.solver.function = 'IPFMINCON';

% opts.method.derivatives = 'internal';
% opts.method.derivatives = 'real-forward';
% opts.method.derivatives = 'real-central';
opts.method.derivatives = 'complex';
% opts.method.derivatives = 'symbolic';

end