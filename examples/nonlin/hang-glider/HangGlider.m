%--------------------------------------------------------------------------
% HangGlider.m
% R. Bulirsch, E. Nerz, H. J. Pesch, and O. von Stryk, "Combining Direct
% and Indirect Methods in Optimal Control: Range Maximization of a Hang
% Glider,*" in Optimal Control, Birkhäuser Basel, 1993, pp. 273–288
% [Online]. Available: http://dx.doi.org/10.1007/978-3-0348-7539-4_20
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = HangGlider(varargin)
% input arguments can be provided in the format 'HangGlider(p,opts)'

% set local functions
ex_opts = @HangGlider_opts; % options function
ex_output = @HangGlider_output; % output function
ex_plot = @HangGlider_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
rho = 1.13;
CD0 = 0.034;
k = 0.069662;
g = 9.80665;
m = 100;
S = 14;
uM = 2.5;
R = 100;
umax = 1.4;
y0 = [0,1000,13.227567500,-1.2875005200];
yf = [inf,900,13.227567500,-1.2875005200];

%% setup
% time horizon
p.t0 = 0; p.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 4; n.np = 1;

% system dynamics
strD{1} = '[';
strD{end+1} = 'p1*y3;';
strD{end+1} = 'p1*y4;';
strD{end+1} = 'p1*( - (S*rho*y3*(k*u1^2 + CD0)*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) - (S*rho*u1*(y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) );';
strD{end+1} = 'p1*( (S*rho*u1*y3*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) - g - (S*rho*(k*u1^2 + CD0)*(y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))*((y4 + uM*exp(-(y1/R - 5/2)^2)*((y1/R - 5/2)^2 - 1))^2 + y3^2)^(1/2))/(2*m) )';
strD{end+1} = ']';
symb.D = horzcat(strD{:});
symb.paramstr = 'rho CD0 k g m S uM R';
symb.param = [rho CD0 k g m S uM R];

% Mayer term
M(1).right = 5; % final states
M(1).left = 0; % singleton
M(1).matrix = [-1,0,0,0];

% simple bounds
UB(1).right = 4; UB(1).matrix = y0'; % initial states
LB(1).right = 4; LB(1).matrix = y0';
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;
UB(3).right = 5; UB(3).matrix = [inf,yf(2:4)]'; % final states
LB(3).right = 5; LB(3).matrix = [-inf,yf(2:4)]';
UB(4).right = 2; UB(4).matrix = [3000 3000 15 15]; % states
LB(4).right = 2; LB(4).matrix = [-3000 0 -15 -15];
UB(5).right = 3; UB(5).matrix = 200; % parameters
LB(5).right = 3; LB(5).matrix = 0;

% guess
Y0 = [[y0];[1250,yf(2:4)]];
U0 = [[1];[1]];
P0 = [[100];[100]];
p.guess = [U0,Y0,P0];

% combine structures
setup.symb = symb; setup.M = M; setup.UB = UB; setup.LB = LB;
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
function opts = HangGlider_opts
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
    opts.qlin.sqpflag = false;
    opts.qlin.trustregionflag = true;
    opts.qlin.improveX0flag = false;
    opts.qlin.delta = 3000;
    opts.qp.disp = 'iter'; % iterations
    opts.qp.solver = 'ipfmincon';
    opts.qlin.method = 'ipfmincon';
    opts.qp.maxiters = 20000;
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 200; % number of nodes
    opts.qlin.sqpflag = false;
    opts.qlin.trustregionflag = true;
    opts.qlin.improveX0flag = false;
    opts.qlin.delta = 3000;
    opts.qp.disp = 'iter'; % iterations
    opts.qp.solver = 'ipfmincon';
    opts.qlin.method = 'ipfmincon';
    opts.qp.maxiters = 20000;
case 3
    opts.general.displevel = 1;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000; % number of nodes
    opts.qlin.method = 'qlin';
    opts.qp.disp = 'none';
    opts.qlin.trustregionflag = false;
    opts.qlin.improveX0flag = false;
end

end