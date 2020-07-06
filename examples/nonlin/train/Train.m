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
% input arguments can be provided in the format 'Train(p,opts)'

% set local functions
ex_opts = @Train_opts;
ex_output = @Train_output;
ex_plot = @Train_plot;

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 4.8;
a = 0.3; b = 0.14; c = 0.16;
z1 = 2; z2 = 4;
s1 = 2; s2 = 0; s3 = -2;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 2; n.ny = 2;

% system dynamics
str{1} = '[';
str{end+1} = 'y2;';
str{end+1} = '(s2-s1)/pi*atan((y1-2)/0.05)+(s3-s2)/pi*atan((y1-4)/0.05)-(a+b*y2+c*y2^2)+u1-u2';
str{end+1} = ']';
symb.D = horzcat(str{:});
symb.paramstr = 'a b c z1 z2 s1 s2 s3';
symb.param = [a b c z1 z2 s1 s2 s3];

% Lagrange term
% symb.Ob = 'u1*y2 + 10^-3*(u1^2 + u2^2)';
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
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.UB = UB; setup.LB = LB;
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