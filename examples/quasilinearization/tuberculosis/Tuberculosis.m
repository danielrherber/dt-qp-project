%--------------------------------------------------------------------------
% Tuberculosis.m
% Sec. 6.16 in J. T. Betts, Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming. SIAM, 2010, isbn: 9780898716887
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Tuberculosis(varargin)
% input arguments can be provided in the format 'Tuberculosis(p,opts)'

% set local functions
ex_opts = @Tuberculosis_opts; % options function
ex_output = @Tuberculosis_output; % output function
ex_plot = @Tuberculosis_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
beta1 = 13;
beta2 = 13;
xmu = 0.0143;
d1 = 0;
d2 = 0;
k1 = 0.5;
k2 = 1;
r1 = 2;
r2 = 1;
p = 0.4;
q = 0.1;
N = 30000;
betas = 0.029;
B1 = 50;
B2 = 500;
lam = xmu*N;

S0 = 76*N/120;
T0 = N/120;
L10 = 36*N/120;
L20 = 2*N/120;
I10 = 4*N/120;
I20 = 1*N/120;

u1min = 0.05;
u1max = 0.95;
u2min = 0.05;
u2max = 0.95;

Nmin = 0;
Nmax = 30000;

%% setup
% time horizon
t0 = 0; tf = 5;

% system dynamics
% y1 = S, y2 = L1, y3 = I1, y4 = L2, y5 = I2, y6 = T
str = {};
str{end+1} = '[';
str{end+1} = 'lam - beta1*y1*y3/N - betas*y1*y5/N - xmu*y1;';
str{end+1} = 'beta1*y1*y3/N - (xmu*k1)*y2 - u1*r1*y2 + (1-u2)*p*r2*y3 + beta2*y6*y3/N - betas*y2*y5/N;';
str{end+1} = 'k1*y2 - (xmu+d1)*y3 - r2*y3;';
str{end+1} = '(1-u2)*q*r2*y3 - (xmu+k2)*y4 + betas*(y1 + y2 + y6)*y5/N;';
str{end+1} = 'k2*y4 - (xmu+d2)*y5;';
str{end+1} = 'u1*r1*y2 + (1 - (1 - u2)*(p+q))*r2*y3 - beta2*y6*y3/N - betas*y6*y5/N - xmu*y6;';
str{end+1} = ']';
str = horzcat(str{:});
symb.D = str;
symb.o.ny = 6; % number of states
symb.o.nu = 2; % number of inputs
symb.o.output = 2; % interp1 compatible
symb.o.param = 'lam beta1 beta2 betas N xmu k1 k2 r1 r2 p q d1 d2';
symb.param = [lam beta1 beta2 betas N xmu k1 k2 r1 r2 p q d1 d2];

% Lagrange term
L(1).left = 1; L(1).right = 1; L(1).matrix = diag(1/2*[B1,B2]); % controls
L(2).left = 0; L(2).right = 2; L(2).matrix = [0,0,0,1,1,0]; % states

% simple bounds
UB(1).right = 4; UB(1).matrix = [S0;L10;I10;L20;I20;T0]; % initial states
LB(1).right = 4; LB(1).matrix = [S0;L10;I10;L20;I20;T0];
UB(2).right = 1; UB(2).matrix = [u1max;u2max]; % controls
LB(2).right = 1; LB(2).matrix = [u1min;u2min];
UB(3).right = 2; UB(3).matrix = repelem(Nmax,6,1); % states
LB(3).right = 2; LB(3).matrix = repelem(Nmin,6,1);

% linear equality constraint
Y(1).linear(1).right = 2; % states
Y(1).linear(1).matrix = [1 1 1 1 1 1]; % states
Y(1).b = N;

% combine structures
setup.symb = symb; setup.L = L; setup.UB = UB; setup.LB = LB; setup.Y = Y;
setup.t0 = t0; setup.tf = tf; setup.p = [];

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
function opts = Tuberculosis_opts
% test number
num = 1;

switch num
case 1
    % default parameters
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 50; % number of nodes
    opts.qlin.sqpflag = false;
    opts.qlin.trustregionflag = true;
    opts.qlin.delta = 2500;
    opts.qp.disp = 'none';
end

end