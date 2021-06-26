%--------------------------------------------------------------------------
% Train.m
% R. J. Vanderbei, "Case Studies in Trajectory Optimization: Trains,
% Planes, and Other Pastimes," Optimization and Engineering, vol. 2, no. 2,
% pp. 215-243, 2001, doi: 10.1023/a:1013145328012
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Train(varargin)
% input arguments can be provided in the format 'Train(auxdata,opts)'

% set local functions
ex_opts = @Train_opts;
ex_output = @Train_output;
ex_plot = @Train_plot;

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 4.8;
a = 0.3; b = 0.14; c = 0.16;
z1 = 2; z2 = 4;
s1 = 2; s2 = 0; s3 = -2;

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = tf;

% number of controls, states, and parameters
n.nu = 2; n.ny = 2;

% system dynamics
str{1} = '[';
str{end+1} = 'y2;';
str{end+1} = '(s2-s1)/pi*atan((y1-2)/0.05)+(s3-s2)/pi*atan((y1-4)/0.05)-(a+b*y2+c*y2^2)+u1-u2';
str{end+1} = ']';
element.dynamics = horzcat(str{:});
element.parameter_list = 'a b c z1 z2 s1 s2 s3';
element.parameter_values = [a b c z1 z2 s1 s2 s3];

% Lagrange term
% element.lagrange = 'u1*y2 + 10^-3*(u1^2 + u2^2)';
L(1).left = 1; L(1).right = 1; L(1).matrix = diag([1e-3,1e-3]);
L(2).left = 1; L(2).right = 2; L(2).matrix = [0,0;1,0;0,0;0,0];
setup.L = L;

% simple bounds
UB(1).right = 4; UB(1).matrix = [0;0]; % initial states
LB(1).right = 4; LB(1).matrix = [0;0];
UB(2).right = 5; UB(2).matrix = [6;0]; % final states
LB(2).right = 5; LB(2).matrix = [6;0];
UB(3).right = 1; UB(3).matrix = [10;2]; % controls
LB(3).right = 1; LB(3).matrix = [0;0];

% guess
Y0 = [[0,0];[6,0]];
U0 = [[10,2];[0,0]];
setup.guess.X = [U0,Y0];

% scaling
setup.scaling(1).right = 1; % controls
setup.scaling(1).matrix = [10,2];
setup.scaling(2).right = 2; % states
setup.scaling(2).matrix = [6,6];

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

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
function opts = Train_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 500; % number of nodes
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
end

end