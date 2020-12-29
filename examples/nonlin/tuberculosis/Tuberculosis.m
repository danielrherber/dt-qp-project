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
xp = 0.4;
q = 0.1;
Npop = 30000;
betas = 0.029;
B1 = 50;
B2 = 500;
lam = xmu*Npop;

S0 = 76*Npop/120;
T0 = Npop/120;
L10 = 36*Npop/120;
L20 = 2*Npop/120;
I10 = 4*Npop/120;
I20 = Npop/120;

u1min = 0.05;
u1max = 0.95;
u2min = 0.05;
u2max = 0.95;

Nmin = 0;
Nmax = Npop;

%% setup
% time horizon
t0 = 0; tf = 5;

% number of controls, states, and parameters
n.nu = 2; n.ny = 6;

% system dynamics
% y1 = S, y2 = L1, y3 = I1, y4 = L2, y5 = I2, y6 = T
str = {};
str{end+1} = '[';
str{end+1} = 'lam - beta1*y1*y3/Npop - betas*y1*y5/Npop - xmu*y1;';
str{end+1} = 'beta1*y1*y3/Npop - (xmu+k1)*y2 - u1*r1*y2 + (1-u2)*xp*r2*y3 + beta2*y6*y3/Npop - betas*y2*y5/Npop;';
str{end+1} = 'k1*y2 - (xmu+d1)*y3 - r2*y3;';
str{end+1} = '(1-u2)*q*r2*y3 - (xmu+k2)*y4 + betas*(y1 + y2 + y6)*y5/Npop;';
str{end+1} = 'k2*y4 - (xmu+d2)*y5;';
str{end+1} = 'u1*r1*y2 + (1 - (1 - u2)*(xp+q))*r2*y3 - beta2*y6*y3/Npop - betas*y6*y5/Npop - xmu*y6;';
str{end+1} = ']';
str = horzcat(str{:});
symb.D = str;
symb.paramstr = 'lam beta1 beta2 betas Npop xmu k1 k2 r1 r2 xp q d1 d2';
symb.param =    [lam beta1 beta2 betas Npop xmu k1 k2 r1 r2 xp q d1 d2];

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
Y(1).b = Npop;

% guess
Y0 = [[S0,L10,I10,L20,I20,T0];[S0,L10,I10,L20,I20,T0]];
U0 = [[0.95,0.95];[0.95,0.95]];
setup.guess.X = [U0,Y0];

% combine structures
setup.symb = symb; setup.L = L; setup.UB = UB; setup.LB = LB; setup.Y = Y;
setup.t0 = t0; setup.tf = tf; setup.p = p; setup.n = n;

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
function opts = Tuberculosis_opts
% test number
num = 2;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.maxiters = 5000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 45; % number of nodes
    opts.solver.tolerance = 1e-12;
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.solver.display = 'iter';
    opts.solver.maxiters = 2000;
case 3 % qlin method
    opts.general.plotflag = 1; % create the plots
    opts.general.displevel = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 200;
    opts.solver.display = 'none';
    opts.method.form = 'qlin';
    opts.method.trustregionflag = true;
    opts.method.sqpflag = false;
    opts.method.delta = inf;
    opts.method.improveguess = false; % disabled
end

end