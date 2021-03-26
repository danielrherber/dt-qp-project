%--------------------------------------------------------------------------
% MineExtraction.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = MineExtraction(varargin)
% input arguments can be provided in the format 'MineExtraction(p,opts)'

% set local functions
ex_opts = @MineExtraction_opts;
ex_output = @MineExtraction_output;
ex_plot = @MineExtraction_plot;

% set p and opts
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 2; % contract length
p.a = 1; % profit rate
p.x0 = 10; % initial ore available
obj_approach = 'string'; % 'function' or 'string'

%% setup
% time horizon
p.t0 = 0; p.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 1;

% system dynamics
setup.A = 0;
setup.B = -1;

% Lagrange term
switch obj_approach
    %----------------------------------------------------------------------
    case 'function'
        symb.Ob = 'u1^2/y1 - a*u1';
    %----------------------------------------------------------------------
    case 'string'
        % provide function, rather than a string, for the objective function
        % NOTE: this feature is currently undocumented and undeveloped
        symb.Ob = []; % only needs to have the field to work
        obj.f = { @(t,param,UYP) UYP(:,1).^2./UYP(:,2) - param(:,1).*UYP(:,1) };
        setup.internalinfo.obj = obj;
    %----------------------------------------------------------------------
end

% problem parameters
symb.paramstr = 'a';
symb.param = [p.a];

% simple bounds
UB(1).right = 4; UB(1).matrix = p.x0;% initial states
LB(1).right = 4; LB(1).matrix = p.x0;
UB(2).right = 2; UB(2).matrix = p.x0;% state bounds
LB(2).right = 2; LB(2).matrix = 0;

% guess
Y0 = [[p.x0];[p.x0]];
U0 = [[0];[0]];
setup.guess.X = [U0,Y0];

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
function opts = MineExtraction_opts
% test number
num = 1;

switch num
    case 1
    % default parameters
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10; % number of nodes
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-10;
    opts.method.form = 'nonlinearprogram';
    opts.method.olqflag = true;
    opts.method.derivatives = 'complex';
end

end