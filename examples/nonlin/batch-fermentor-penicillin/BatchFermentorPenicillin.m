%--------------------------------------------------------------------------
% BatchFermentorPenicillin.m
% J. R. Banga, E. Balsa-Canto, C. G. Moles, and A. A. Alonso, "Dynamic
% optimization of bioprocesses: Efficient and robust numerical strategies,"
% Journal of Biotechnology, vol. 117, no. 4, pp. 407–419, Jun. 2005,
% doi: 10.1016/j.jbiotec.2005.02.013
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BatchFermentorPenicillin(varargin)
% input arguments can be provided in the format 'BatchFermentorPenicillin(p,opts)'

% set local functions
ex_opts = @BatchFermentorPenicillin_opts;
ex_output = @BatchFermentorPenicillin_output;
ex_plot = @BatchFermentorPenicillin_plot;

% set p and opts
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters

%% setup
% time horizon (scaled)
p.t0 = 0; p.tf = 1;

% number of controls, states, and parameters
n.nu = 1; n.ny = 4; n.np = 1;

% system dynamics
str{1} = '[';
str{end+1} = 'p1*( (0.11*(y3./(0.006*y1+y3)).*y1-u1.*(y1./500./y4)) ); ';
str{end+1} = 'p1*( (0.0055*(y3./(0.0001+y3.*(1+10*y3))).*y1-0.01*y2-u1.*(y2./500./y4)) ); ';
str{end+1} = 'p1*( (-0.11*(y3./(0.006*y1+y3)).*y1/0.47-0.0055*(y3./(0.0001+y3.*(1+10*y3))).*y1/1.2-y1.*(0.029*y3./(0.0001+y3))+u1./y4.*(1-y3/500)) ); ';
str{end+1} = 'p1*( u1/500 )';
str{end+1} = ']';
symb.D = horzcat(str{:});

% Mayer term
Mmatrix = zeros(4); Mmatrix(2,4) = -1;
M(1).left = 5; M(1).right = 5; M(1).matrix = Mmatrix;

% simple bounds
UB(1).right = 4; UB(1).matrix = [1.5,0,0,7]; % initial states
LB(1).right = 4; LB(1).matrix = [1.5,0,0,7];
UB(2).right = 1; UB(2).matrix = 50; % controls
LB(2).right = 1; LB(2).matrix = 0;
UB(3).right = 5; UB(3).matrix = [40,50,25,10]; % final states
LB(3).right = 5; LB(3).matrix = [0,0,0,0];
UB(4).right = 2; UB(4).matrix = [40,50,25,10]; % states
LB(4).right = 2; LB(4).matrix = [0,0,0,0];

% guess
Y0 = [[1.5,0,0,7];[30,8.5,0,10]];
U0 = [[10];[10]];
P0 = [[10];[10]];
setup.guess.X = [U0,Y0,P0];

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
function opts = BatchFermentorPenicillin_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 150; % number of nodes
    opts.solver.maxiters = 20000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.solver.tolerance = 1e-8;
    opts.method.derivatives = 'symbolic';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100; % number of nodes
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivativeflag = true;
end

end