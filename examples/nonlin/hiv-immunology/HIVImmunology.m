%--------------------------------------------------------------------------
% HIVImmunology.m
% H. R. Joshi, "Optimal control of an HIV immunology model," Optimal
% Control Applications and Methods, vol. 23, no. 4, pp. 199–213, 2002,
% doi: 10.1002/oca.710
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = HIVImmunology(varargin)
% input arguments can be provided in the format 'HIVImmunology(p,opts)'

% set local functions
ex_opts = @HIVImmunology_opts;
ex_output = @HIVImmunology_output;
ex_plot = @HIVImmunology_plot;

% set p and opts
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 50;
s1 = 2; s2 = 1.5;
m = 0.002; C = 0.007;
k = 2.5e-4; g = 30;
b1 = 14; b2 = 1;
A1 = 2.5e5; A2 = 75;

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 2; n.ny = 2;

% system dynamics
str{1} = '[';
str{end+1} = 's1-s2*y2/(b1+y2)- m*y1-k*y1*y2+u1*y1; ';
str{end+1} = 'g*(1-u2)*y2/(b2+y2)-C*y1*y2';
str{end+1} = ']';
symb.D = horzcat(str{:});
symb.paramstr = 's1 s2 m C k g b1 b2 A1 A2';
symb.param = [s1 s2 m C k g b1 b2 A1 A2];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = diag([A1 A2]);
L(2).left = 0; L(2).right = 2; L(2).matrix = -[1 0];

% simple bounds
UB(1).right = 4; UB(1).matrix = [400;3];% initial states
LB(1).right = 4; LB(1).matrix = [400;3];
UB(2).right = 2; UB(2).matrix = [1200;5]; % states
LB(2).right = 2; LB(2).matrix = [0;0.05];
UB(3).right = 1; UB(3).matrix = [0.02;0.9]; % controls
LB(3).right = 1; LB(3).matrix = [0;0];

% guess
Y0 = [[400,3];[400,3]];
U0 = [[0.02,0.9];[0,0]];
p.guess = [U0,Y0];

% combine structures
setup.symb = symb; setup.UB = UB; setup.LB = LB; setup.L = L;
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
function opts = HIVImmunology_opts
% test number
num = 1;

switch num
    case 1
    % default parameters
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
end

end