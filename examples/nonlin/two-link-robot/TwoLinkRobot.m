%--------------------------------------------------------------------------
% TwoLinkRobot.m
% Section 12.4.2 of R. Luus, Iterative Dynamic Programming. CRC Press,
% 2000, isbn: 1584881488
%--------------------------------------------------------------------------
% Similar to http://www.ee.ic.ac.uk/ICLOCS/ExampleRobotArm.html, but the
% equations at this link have some differences from Luus2000
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = TwoLinkRobot(varargin)
% input arguments can be provided in the format 'TwoLinkRobot(p,opts)'

% set local functions
ex_opts = @TwoLinkRobot_opts; % options function
ex_output = @TwoLinkRobot_output; % output function
ex_plot = @TwoLinkRobot_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
y0 = [0 0 0.5 0]; % initial states
yf = [0 0 0.5 0.522]; % final states

%% setup
% time horizon (scaled)
p.t0 = 0; p.tf = 1;

% number of controls, states, and parameters
n.nu = 2; n.ny = 4; n.np = 1;

% system dynamics
str{1} = '[';
str{2} = 'p1*( (sin(y3).*(9/4*cos(y3).*y1.^2+2*y2.^2) + 4/3*(u1-u2) - 3/2*cos(y3).*u2 )./ (31/36 + 9/4*sin(y3).^2) ); ';
str{3} = 'p1*( -( sin(y3).*(9/4*cos(y3).*y2.^2+7/2*y1.^2) - 7/3*u2 + 3/2*cos(y3).*(u1-u2) )./ (31/36 + 9/4*sin(y3).^2) ); ';
str{4} = 'p1*( y2-y1 ); ';
str{5} = 'p1*( y1 )';
str{6} = ']';
element.dynamics = horzcat(str{:});

% Mayer term
M(1).right = 3; M(1).left = 0; M(1).matrix = 1; % parameters

% Lagrange term
% L(1).right = 1; L(1).left = 1; L(1).matrix = 0.01*eye(2); % controls

% simple bounds
UB(1).right = 4; UB(1).matrix = y0; % initial states
LB(1).right = 4; LB(1).matrix = y0;
UB(2).right = 5; UB(2).matrix = yf; % final states
LB(2).right = 5; LB(2).matrix = yf;
UB(3).right = 1; UB(3).matrix = [1 1]; % controls
LB(3).right = 1; LB(3).matrix = [-1 -1];
LB(4).right = 3; LB(4).matrix = 0.1; % parameters

% guess
Y0 = [y0;yf];
U0 = [[1 1];[-1 -1]];
P0 = [[3.1];[3.1]];
setup.guess.X = [U0,Y0,P0];

% combine structures
setup.element = element; setup.M = M;  setup.UB = UB; setup.LB = LB; % setup.L = L;
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
function opts = TwoLinkRobot_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 200; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
end

end