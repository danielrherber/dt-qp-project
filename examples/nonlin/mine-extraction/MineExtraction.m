%--------------------------------------------------------------------------
% MineExtraction.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = MineExtraction(varargin)
% input arguments can be provided in the format 'MineExtraction(auxdata,opts)'

% set local functions
ex_opts = @MineExtraction_opts;
ex_output = @MineExtraction_output;
ex_plot = @MineExtraction_plot;

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
tf = 2; % contract length
auxdata.a = 1; % profit rate
auxdata.x0 = 10; % initial ore available
obj_approach = 'function'; % 'function' or 'string'

%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = tf;

% number of controls, states, and parameters
n.nu = 1; n.ny = 1;

% system dynamics
setup.A = 0;
setup.B = -1;

% Lagrange term
switch obj_approach
    %----------------------------------------------------------------------
    case 'string'
        element.lagrange = 'u1^2/y1 - a*u1';
    %----------------------------------------------------------------------
    case 'function'
        % provide function, rather than a string, for the objective function
        % NOTE: this feature is currently undocumented and undeveloped
        element.lagrange = []; % only needs to have the field to work
        obj.f = { @(t,param,UYP) UYP(:,1).^2./UYP(:,2) - param(:,1).*UYP(:,1) };
        setup.internalinfo.obj = obj;
    %----------------------------------------------------------------------
end

% problem parameters
element.parameter_list = 'a';
element.parameter_values = [auxdata.a];

% simple bounds
UB(1).right = 4; UB(1).matrix = auxdata.x0;% initial states
LB(1).right = 4; LB(1).matrix = auxdata.x0;
UB(2).right = 2; UB(2).matrix = auxdata.x0;% state bounds
LB(2).right = 2; LB(2).matrix = 0;

% guess
Y0 = [[auxdata.x0];[auxdata.x0]];
U0 = [[0];[0]];
setup.guess.X = [U0,Y0];

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