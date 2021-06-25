%--------------------------------------------------------------------------
% GreenhouseClimate.m
% pp. 15-24 of [1]G. van Straten, G. van Willigenburg, E. van Henten, and
% R. van Ooteghem, Optimal Control of Greenhouse Cultivation. CRC Press,
% 2010 [Online]. Available: http://dx.doi.org/10.1201/b10321
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = GreenhouseClimate(varargin)
% input arguments can be provided in the format 'GreenhouseClimate(p,opts)'

% set local functions
ex_opts = @GreenhouseClimate_opts; % options function
ex_output = @GreenhouseClimate_output; % output function
ex_plot = @GreenhouseClimate_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% setup
casenum = 1;

switch casenum
    %----------------------------------------------------------------------
    case 1 % textbook example: LQDO
    % problem parameters
    xtf = 48;
    xp1 = 7.5e-8;
    xp2 = 1;
    xp3 = 0.1;
    xp4 = 4.55e-4;
    xp5 = 136.4;
    umax = 100;

    % system dynamics
    A = cell(2);
    A{1,1} = 0;
    A{1,2} = @(t) xp1*( max(0,800*sin(4*pi*t/xtf - 0.65*pi)) );
    A{2,1} = 0;
    A{2,2} = -xp2;
    setup.A = A;
    setup.B = [0;xp3];
    d = cell(2,1);
    d{2,1} = @(t) xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi));
    setup.d = d;

    opts.method.form = 'lqdo';
    %----------------------------------------------------------------------
    case 2 % PROPT example: LQDO
    % problem parameters
    xtf = 48;
    xp1 = 3e-6/40;
    xp2 = 1;
    xp3 = 0.1;
    xp4 = 7.5e-2/220;
    xp5 = 3e4/220;
    umax = 10;

    % system dynamics
    A = cell(2);
    A{1,1} = 0;
    A{1,2} = @(t) xp1*( 800*sin(4*pi*t/xtf - 0.65*pi) );
    A{2,1} = 0;
    A{2,2} = -xp2;
    setup.A = A;
    setup.B = [0;xp3];
    d = cell(2,1);
    d{2,1} = @(t) xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi));
    setup.d = d;

    opts.method.form = 'lqdo';
    %----------------------------------------------------------------------
    case 3 % textbook example
    % problem parameters
    xtf = 48;
    xp1 = 7.5e-8;
    xp2 = 1;
    xp3 = 0.1;
    xp4 = 4.55e-4;
    xp5 = 136.4;
    umax = 100;

    % system dynamics
    str{1} = '[';
    str{end+1} = '1/2*xp1*( (800*sin(4*pi*t/xtf - 0.65*pi) + abs(800*sin(4*pi*t/xtf - 0.65*pi)) ) )*y2; '; % max function
    str{end+1} = 'xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi) - y2) + xp3*u1';
    str{end+1} = ']';
    element.dynamics = horzcat(str{:});
    element.parameter_list = 'xp1 xp2 xp3 xp4 xp5 xtf';
    element.parameter_values = [xp1 xp2 xp3 xp4 xp5 xtf];
    setup.element = element;

    opts.method.form = 'nonlinearprogram';
    %----------------------------------------------------------------------
    case 4 % PROPT example
    % problem parameters
    xtf = 48;
    xp1 = 3e-6/40;
    xp2 = 1;
    xp3 = 0.1;
    xp4 = 7.5e-2/220;
    xp5 = 3e4/220;
    umax = 10;

    % system dynamics
    str{1} = '[';
    str{end+1} = 'xp1*( 800*sin(4*pi*t/xtf - 0.65*pi) )*y2; ';
    str{end+1} = 'xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi) - y2) + xp3*u1';
    str{end+1} = ']';
    element.dynamics = horzcat(str{:});
    element.parameter_list = 'xp1 xp2 xp3 xp4 xp5 xtf';
    element.parameter_values = [xp1 xp2 xp3 xp4 xp5 xtf];
    setup.element = element;

    opts.method.form = 'nonlinearprogram';
end

p.casenum = casenum; p.xtf = xtf; p.xp1 = xp1; p.xp2 = xp2; p.xp3 = xp3;
p.xp4 = xp4; p.xp1 = xp5;

% number of controls, states, and parameters
n.nu = 1; n.ny = 2;

% time horizon
p.t0 = 0; p.tf = xtf;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-xp5,0];

% Lagrange term
L(1).left = 0; L(1).right = 1; L(1).matrix = xp4;

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,10]; % initial states
LB(1).right = 4; LB(1).matrix = [0,10];
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;

% combine structures
setup.M = M; setup.L = L; setup.UB = UB; setup.LB = LB;
setup.t0 = p.t0; setup.tf = p.tf; setup.p = p; setup.n = n;

%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
ex_plot(T,U,Y,P,F,in,opts,sol)

end
% User options function for this example
function opts = GreenhouseClimate_opts
% test number
num = 1;

switch num
case 1
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 241; % number of nodes
    opts.solver.tolerance = 1e-8;
    opts.solver.maxiters = 200;
    opts.solver.display = 'iter';
end

end