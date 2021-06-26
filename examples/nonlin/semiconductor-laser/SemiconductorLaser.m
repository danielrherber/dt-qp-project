%--------------------------------------------------------------------------
% SemiconductorLaser.m
% J.-H. R. Kim, G. L. Lippi, and H. Maurer, "Minimizing the transition time
% in lasers by optimal control methods," Physica D: Nonlinear Phenomena,
% vol. 191, no. 3–4, pp. 238–260, May 2004, doi: 10.1016/j.physd.2003.12.002
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = SemiconductorLaser(varargin)
% input arguments can be provided in the format 'SemiconductorLaser(auxdata,opts)'

% set local functions
ex_opts = @SemiconductorLaser_opts;
ex_output = @SemiconductorLaser_output;
ex_plot = @SemiconductorLaser_plot;

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
taup = 2.072e-12;
gam = 0.3;
Gp = 2.628e4;
Ntr = 7.8e7;
q = 1.60219e-19;
ep = 9.6e-8;
beta = 1.735e-4;
A = 1e8;
B = 2.788;
C = 7.3e-9;
P0 = 1.5e7;
S0 = 0.6119512914e5;
Sf = 3.4063069073e5;
N0 = 1.3955581328e8;
Nf = 1.4128116637e8;
Imin = 2e-3;
Imax = 67.5e-3;
ts = 1e-12;

% assign to parameter structure
auxdata.S0 = S0;
auxdata.N0 = N0;
auxdata.Imin = Imin;

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2; n.np = 1;

% system dynamics
str = {};
str{end+1} = '[';
str{end+1} = 'ts*p1/S0*( -y1*S0/taup + gam*Gp*(y2*N0-Ntr)*(1-ep*y1*S0)*y1*S0 + beta_*B_*y2*N0*(y2*N0 + P0));';
str{end+1} = 'ts*p1/N0*( Imin*u1/q - (A_*y2*N0 + B_*y2*N0*(y2*N0+P0) + C_*y2*N0*(y2*N0+P0)^2) - gam*Gp*(y2*N0-Ntr)*(1-ep*y1*S0)*y1*S0);';
str{end+1} = ']';
str = horzcat(str{:});
element.dynamics = str;

% problem parameters
element.parameter_list = 'taup ep P0 B_ Imin Gp Ntr q C_ gam beta_ A_ S0 N0 ts';
element.parameter_values = [taup ep P0 B Imin Gp Ntr q C gam beta A S0 N0 ts];

% Mayer term
M(1).left = 0; M(1).right = 3; M(1).matrix = [1];

% simple bounds
UB(1).right = 4; UB(1).matrix = [S0/S0,N0/N0]; % initial states
LB(1).right = 4; LB(1).matrix = [S0/S0,N0/N0];
UB(2).right = 5; UB(2).matrix = [Sf/S0,Nf/N0]; % final states
LB(2).right = 5; LB(2).matrix = [Sf/S0,Nf/N0];
UB(3).right = 1; UB(3).matrix = Imax/Imin; % controls
LB(3).right = 1; LB(3).matrix = Imin/Imin;
UB(4).right = 3; UB(4).matrix = 100; % parameters
LB(4).right = 3; LB(4).matrix = 0.1;

% guess
Y0 = [[S0/S0,N0/N0];[Sf/S0,Nf/N0]];
U0 = [[Imin/Imin];[Imax/Imin]];
P0 = [[100];[100]];
setup.guess.X = [U0,Y0,P0];

% combine structures
setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
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
function opts = SemiconductorLaser_opts
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
    opts.dt.nt = 1000; % number of nodes
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
end

end