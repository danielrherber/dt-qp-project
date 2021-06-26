%--------------------------------------------------------------------------
% Multiextremal.m
% Problems 4-12 in A. Yu. Gornov, T. S. Zarodnyuk, T. I. Madzhara, A. V.
% Daneeva, and I. A. Veyalko, "A Collection of Test Multiextremal Optimal
% Control Problems*", in Optimization, Simulation, and Control, Springer
% New York, 2012, pp. 257–274, doi: 10.1007/978-1-4614-5131-0_16
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Multiextremal(varargin)
% input arguments can be provided in the format 'Multiextremal(auxdata,opts)'

% set local functions
ex_opts = @Multiextremal_opts; % options function
ex_output = @Multiextremal_output; % output function
ex_plot = @Multiextremal_plot; % plot function

% set auxdata and opts (see local_opts)
[auxdata,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% tunable parameters
testnum = 6; % 4-12 from the reference above
switch testnum
    %----------------------------------------------------------------------
    case 4 % global solution found
    Dy1 = 'y2';
    Dy2 = 'u1 - sin(y1)';
    tf = 5;
    y0 = [5,0];
    umin = -1; umax = 1;

    % Mayer term
    M(1).left = 5; M(1).right = 5; M(1).matrix = diag([1,1]); % final states^2
    %----------------------------------------------------------------------
    case 5 % global solution found
    Dy1 = 'sin(y2)';
    Dy2 = 'u1 - exp(y1)';
    tf = 5;
    y0 = [1,1];
    umin = -1; umax = 1;

    % Mayer term
    M(1).left = 0; M(1).right = 5; M(1).matrix = [1,1]; % final states
    %----------------------------------------------------------------------
    case 6 % global solution not found
    Dy1 = '1 - y2^2 + 0.5*y2';
    Dy2 = 'y1*u1';
    tf = 1.7;
    y0 = [0.5,-0.2];
    umin = -2; umax = 2;

    % Mayer term
    M(1).left = 5; M(1).right = 5; M(1).matrix = [0,1;0,0]; % final states^2
    M(2).left = 0; M(2).right = 5; M(2).matrix = [1,-1]; % final states
    %----------------------------------------------------------------------
    case 7 % global solution found
    Dy1 = 'sin(y2) + cos(t)';
    Dy2 = '(u1 + t)^2';
    tf = 1.5;
    y0 = [1,1];
    umin = -3; umax = 3;

    % Mayer term
    M(1).left = 0; M(1).right = 5; M(1).matrix = [1,1]; % final states
    %----------------------------------------------------------------------
    case 8 % global solution found
    Dy1 = '-(2+u1)*(y1+0.25) + (y2+0.5)*exp(25*y1/(y1+2))';
    Dy2 = '0.5 - y2 - (y2+0.5)*exp(25*y1/(y1+2))';
    tf = 0.78;
    y0 = [0.09,0.09];
    umin = 0; umax = 5;

    % Lagrange term
    L(1).left = 1; L(1).right = 1; L(1).matrix = 0.1; % controls^2
    L(2).left = 2; L(2).right = 2; L(2).matrix = diag([1,1]); % states^2
    %----------------------------------------------------------------------
    % case 9 % nonlinear Mayer term
    %----------------------------------------------------------------------
    % case 10 % nonlinear Mayer term
    %----------------------------------------------------------------------
    case 11 % global solution found
    Dy1 = 'cos(y2)';
    Dy2 = 'u1 - sin(y1)';
    tf = 7;
    y0 = [0,0];
    umin = -1; umax = 1;

    % Mayer term
    M(1).left = 0; M(1).right = 5; M(1).matrix = [1,0]; % final states
    %----------------------------------------------------------------------
    case 12 % global solution not found
    Dy1 = 'y2';
    Dy2 = 'u1 - y1 + y1^3/6 - y1^5/120';
    tf = 7;
    y0 = [5,0];
    umin = -1; umax = 1;

    % Mayer term
    M(1).left = 5; M(1).right = 5; M(1).matrix = eye(2); % final states^2
    M(2).left = 0; M(2).right = 5; M(2).matrix = [2,-19]; % final states
    M(3).left = 0; M(3).right = 0; M(3).matrix = 365/4; % constant
    %----------------------------------------------------------------------
end
%% setup
% time horizon
auxdata.t0 = 0; auxdata.tf = tf;

% number of controls, states, and parameters
n.ny = 2; n.nu = 1;

% system dynamics
str{1} = '[';
str{end+1} = Dy1;
str{end+1} = '; ';
str{end+1} = Dy2;
str{end+1} = ']';
element.dynamics = horzcat(str{:});

% simple bounds
UB(1).right = 4; UB(1).matrix = y0; % initial states
LB(1).right = 4; LB(1).matrix = y0;
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = umin;

% guess
Y0 = [[y0];[y0]];
U0 = [[0];[0]];
setup.guess.X = [U0,Y0];

% combine structures
setup.element = element; setup.UB = UB; setup.LB = LB;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;
if exist('M','var')
	setup.M = M;
end
if exist('L','var')
	setup.L = L;
end

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
function opts = Multiextremal_opts
% test number
num = 1;

switch num
case 1
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 1000; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-12;
    opts.method.form = 'nonlinearprogram';
case 2
    opts.general.displevel = 2;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100; % number of nodes
    opts.solver.display = 'iter'; % iterations
    opts.solver.function = 'ipfmincon';
    opts.solver.tolerance = 1e-12;
    opts.method.form = 'nonlinearprogram';
    opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
    opts.dt.meshr.tolerance = 1e-7;
end

end