%--------------------------------------------------------------------------
% TransferMinFuel.m
% I. M. Ross, Q. Gong, and P. Sekhavat, "Low-Thrust, High-Accuracy
% Trajectory Optimization," Journal of Guidance, Control, and Dynamics,
% vol. 30, no. 4, pp. 921–933, Jul. 2007, doi: 10.2514/1.23181
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clearvars -except externalInput

% set options and standardize
if ~exist('externalInput','var')
    opts = localOpts;
end
DTQP_standardizedinputs2

% create setup structure
setup = createSetup;

% solve with DTQP
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% post-processing
[O,sol] = postProcessing(T,U,Y,P,F,in,opts);

%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% auxiliary data
testnum = 22; % see below
r0 = 1; rf = 4; % initial and final radius
auxdata.tf = 57; % final time
umax = 0.01; % maximum thrust
setup.auxdata = auxdata;

% scaling
scalenum = 2;
switch scalenum
    %----------------------------------------------------------------------
    case 1
    Us = 1; Ts = 1; Rs = 1; As = 1; Vrs = 1;
    %----------------------------------------------------------------------
    case 2
	Us = umax; Ts = auxdata.tf; Rs = (r0+rf)/2; As = 10; Vrs = 0.1;
end
auxdata.Us = Us; auxdata.Ts = Ts; auxdata.Rs = Rs; auxdata.As = As; auxdata.Vrs = 0.1;

% time horizon
setup.t0 = 0; setup.tf = auxdata.tf/Ts;

% number of controls, states, and parameters
setup.counts.nu = 2;
setup.counts.nx = 4;

% simple bounds
Ys = [Rs As Vrs 1];
auxdata.Ys = Ys;
setup.auxdata = auxdata;
UB(1).right = 4; UB(1).matrix = [r0,0,0,1]./Ys; % initial states
LB(1).right = 4; LB(1).matrix = [r0,0,0,1]./Ys;
UB(2).right = 5; UB(2).matrix = [rf,inf,0,0.5]./Ys; % final states
LB(2).right = 5; LB(2).matrix = [rf,-inf,0,0.5]./Ys;

% system dynamics
str{1} = '[';
str{end+1} = 'Ts*Vrs/Rs*y3; ';
str{end+1} = 'Ts/As/Rs*y4/y1; ';
str{end+1} = 'Ts/Vrs/Rs*y4^2/y1 - Ts/Vrs/Rs^2/y1^2 + Ts/Vrs*Us*u1; ';
str{end+1} = '-Ts*Vrs/Rs*y3*y4/y1 + Ts*Us*u2';
str{end+1} = ']';
setup.nonlin.dynamics = horzcat(str{:});

% symbolic data for nonlin
setup.nonlin.data.symbols = 'umax Us Ts Rs As Vrs';
setup.nonlin.data.values = [umax Us Ts Rs As Vrs];

switch testnum
    %----------------------------------------------------------------------
    case {1,11} % p = 1, q = inf

    if testnum == 1 % original formulation
        % Lagrange term
        % setup.nonlin.lagrange = 'abs(u1) + abs(u2)';
        setup.nonlin.lagrange = 'sqrt(u1^2) + sqrt(u2^2)'; % works with complex numbers

        % linear inequality constraint
        UB(3).right = 1; UB(3).matrix = [umax, umax]/Us; % controls
        LB(3).right = 1; LB(3).matrix = [-umax, -umax]/Us;

        % guess
        Y0 = [[r0,0,0,1];[rf,0,0,0.5]]./Ys;
        U0 = [[0.01,0.01];[0.01,0.01]]/Us;
        setup.method.guess.X = [U0,Y0];

    elseif testnum == 11 % transformed problem
        % two additional controls
        setup.counts.nu = setup.counts.nu + 2;

        % Lagrange term
        % element.lagrange = 'u3 + u4';
        L(1).left = 0; L(1).right = 1; L(1).matrix = [0 0 1 1];
        setup.lq.lagrange = L;

        % general linear inequality constraints
        Z(1).linear(1).right = 1; Z(1).linear(1).matrix = [1 0 -1 0];
        Z(1).b = 0; % u1 - u3 < 0
        Z(2).linear(1).right = 1; Z(2).linear(1).matrix = [-1 0 -1 0];
        Z(2).b = 0; % -u1 - u3 < 0
        Z(3).linear(1).right = 1; Z(3).linear(1).matrix = [0 1 0 -1];
        Z(3).b = 0; % u2 - u4 < 0
        Z(4).linear(1).right = 1; Z(4).linear(1).matrix = [0 -1 0 -1];
        Z(4).b = 0; % -u2 - u4 < 0
        setup.lq.inequality = Z;

        % linear inequality constraint
        UB(3).right = 1; UB(3).matrix = [umax, umax, umax, umax]/Us; % controls
        LB(3).right = 1; LB(3).matrix = [-umax, -umax, 0, 0]/Us;

        % guess
        Y0 = [[r0,0,0,1];[rf,0,0,0.5]]./Ys;
        U0 = [[0,umax,0,umax];[0,umax,0,umax]]/Us;
        setup.method.guess.X = [U0,Y0];

    end

    %----------------------------------------------------------------------
    case {2,22} % p = 2, q = 2

    if testnum == 2 % original formulation
        % Lagrange term
        setup.nonlin.lagrange = '(u1^2 + u2^2)^(1/2)';

        % nonlinear inequality constraint
        % setup.nonlin.inequality.func = '(u1^2 + u2^2)^(1/2) - umax/Us';
        setup.nonlin.inequality.func = '(u1^2 + u2^2) - umax^2/Us^2'; % better form
        setup.nonlin.inequality.pathboundary = 1;

        % linear inequality constraint
        UB(3).right = 1; UB(3).matrix = [umax, umax]/Us; % controls
        LB(3).right = 1; LB(3).matrix = [-umax, -umax]/Us;

        % guess
        Y0 = [[r0,0,0,1];[rf,0,0,0.5]]./Ys;
        U0 = [[0,umax];[0,umax]]/Us;
        setup.method.guess.X = [U0,Y0];

    elseif testnum == 22 % transformed problem
        % one additional control
        setup.counts.nu = setup.counts.nu + 1;

        % Lagrange term
        L(1).left = 0; L(1).right = 1; L(1).matrix = [0,0,1];
        setup.lq.lagrange = L;

        % nonlinear inequality constraint
        % setup.nonlin.equality.func = 'u1^2 + u2^2 - u3^2';
        % setup.nonlin.equality.pathboundary = 1;
        setup.nonlin.inequality.func = 'u1^2 + u2^2 - u3^2'; % better form
        setup.nonlin.inequality.pathboundary = 1;

        % linear inequality constraint
        UB(3).right = 1; UB(3).matrix = [1.05*umax, 1.05*umax, umax]/Us; % controls
        LB(3).right = 1; LB(3).matrix = [-1.05*umax, -1.05*umax, 0]/Us;

        % guess
        Y0 = [[r0,0,0,1];[rf,20,0,0.5]]./Ys;
        U0 = [[0,umax,umax];[0,umax,umax]]/Us;
        setup.method.guess.X = [U0,Y0];

    end
end

% simple bounds
setup.lq.ub = UB; setup.lq.lb = LB;

end

%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,opts)

% unscale
T = T*in.auxdata.Ts;
U = U*in.auxdata.Us;
Y = Y.*repmat(in.auxdata.Ys,in.nt,1);
F = F*in.auxdata.Ts*in.auxdata.Us;

% outputs
[O,sol] = TransferMinFuel_output(T,U,Y,P,F,in,opts);

% plots
TransferMinFuel_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts

% test number
num = 1;

switch num
case 1 % ipfmincon method
    opts.general.plotflag = 1; % create the plots
    opts.general.saveflag = false;
    opts.general.displevel = 2;
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 300; % number of nodes
    opts.solver.tolerance = 1e-10;
    opts.solver.maxiters = 2000;
    opts.solver.display = 'iter';
    opts.solver.function = 'ipfmincon';
    opts.method.form = 'nonlinearprogram';
    opts.olqflag = true;
end

end