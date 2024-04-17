%--------------------------------------------------------------------------
% GreenhouseClimate.m
% pp. 15-24 of [1]G. van Straten, G. van Willigenburg, E. van Henten, and
% R. van Ooteghem, Optimal Control of Greenhouse Cultivation. CRC Press,
% 2010 [Online]. Available: http://dx.doi.org/10.1201/b10321
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
[O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts);


%--------------------------------------------------------------------------
% create setup structure
%--------------------------------------------------------------------------
function setup = createSetup

% initialize struct
setup = DTQP_setup_initialize;

% tunable parameters
casenum = 4;

switch casenum

    case 1 % textbook example: LQDO

        % problem parameters
        xtf = 48;      xp1 = 7.5e-8; xp2 = 1;    xp3 = 0.1;
        xp4 = 4.55e-4; xp5 = 136.4;  umax = 100;

        % system dynamics
        A = cell(2);
        A{1,1} = 0; A{1,2} = @(t) xp1*( max(0,800*sin(4*pi*t/xtf - 0.65*pi)));
        A{2,1} = 0; A{2,2} = -xp2;
        setup.lq.dynamics.A = A;
        setup.lq.dynamics.B = [0;xp3];
        d = cell(2,1);
        d{2,1} = @(t) xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi));
        setup.lq.dynamics.Bv = d;
        % opts.method.form = 'lqdo';

    case 2 % PROPT example: LQDO

        % problem parameters
        xtf = 48;         xp1 = 3e-6/40; xp2 = 1;   xp3 = 0.1;
        xp4 = 7.5e-2/220; xp5 = 3e4/220; umax = 10;

        % system dynamics
        A = cell(2);
        A{1,1} = 0; A{1,2} = @(t) xp1*(800*sin(4*pi*t/xtf - 0.65*pi));
        A{2,1} = 0; A{2,2} = -xp2;
        setup.lq.dynamics.A = A;
        setup.lq.dynamics.B = [0;xp3];
        d = cell(2,1);
        d{2,1} = @(t) xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi));
        setup.lq.dynamics.Bv = d;
        % opts.method.form = 'lqdo';

    case 3 % textbook example

        % problem parameters
        xtf = 48; xp1 = 7.5e-8; xp2 = 1; xp3 = 0.1;
        xp4 = 4.55e-4; xp5 = 136.4; umax = 100;

        % number of controls, states, and parameters
        setup.counts.nu = 1; setup.counts.nx = 2;

        % system dynamics
        str{1} = '[';
        str{end+1} = '1/2*xp1*( (800*sin(4*pi*t/xtf - 0.65*pi) + abs(800*sin(4*pi*t/xtf - 0.65*pi))))*y2; '; % max function
        str{end+1} = 'xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi) - y2) + xp3*u1';
        str{end+1} = ']';
        setup.nonlin.dynamics = horzcat(str{:});
        setup.nonlin.data.symbols = 'xp1 xp2 xp3 xp4 xp5 xtf';
        setup.nonlin.data.values = [xp1 xp2 xp3 xp4 xp5 xtf];
        % opts.method.form = 'nonlinearprogram';

    case 4 % PROPT example

        % problem parameters
        xtf = 48; xp1 = 3e-6/40; xp2 = 1; xp3 = 0.1;
        xp4 = 7.5e-2/220; xp5 = 3e4/220; umax = 10;

        % number of controls, states, and parameters
        setup.counts.nu = 1; setup.counts.nx = 2;


        % system dynamics
        str{1} = '[';
        str{end+1} = 'xp1*( 800*sin(4*pi*t/xtf - 0.65*pi) )*y2; ';
        str{end+1} = 'xp2*( 15 + 10*sin(4*pi*t/xtf - 0.65*pi) - y2) + xp3*u1';
        str{end+1} = ']';
        setup.nonlin.dynamics = horzcat(str{:});
        setup.nonlin.data.symbols = 'xp1 xp2 xp3 xp4 xp5 xtf';
        setup.nonlin.data.values = [xp1 xp2 xp3 xp4 xp5 xtf];
        % opts.method.form = 'nonlinearprogram';
end


auxdata.casenum = casenum; auxdata.xtf = xtf; auxdata.xp1 = xp1;
auxdata.xp2 = xp2; auxdata.xp3 = xp3; auxdata.xp4 = xp4; auxdata.xp1 = xp5;
setup.auxdata = auxdata;

% time horizon
setup.t0 = 0; setup.tf = xtf;

% Lagrange term
L(1).left = 0; L(1).right = 1; L(1).matrix = xp4;
setup.lq.lagrange = L;

% Mayer term
M(1).left = 0; M(1).right = 5; M(1).matrix = [-xp5,0];
setup.lq.mayer = M;

% simple bounds
UB(1).right = 4; UB(1).matrix = [0,10]; % initial states
LB(1).right = 4; LB(1).matrix = [0,10];
UB(2).right = 1; UB(2).matrix = umax; % controls
LB(2).right = 1; LB(2).matrix = 0;
setup.lq.ub = UB; setup.lq.lb = LB;

end



%--------------------------------------------------------------------------
% post-processing
%--------------------------------------------------------------------------
function [O,sol] = postProcessing(T,U,Y,P,F,in,setup,opts)

% outputs
[O,sol] = GreenhouseClimate_output(T,U,Y,P,F,in,opts);

% plots
GreenhouseClimate_plot(T,U,Y,P,F,in,opts,sol)

end

%--------------------------------------------------------------------------
% local options
%--------------------------------------------------------------------------
function opts = localOpts
% test number
num = 1;

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

% Case 1 and 2
% opts.method.form = 'lqdo';

% case 3 & 4
% opts.method.form = 'nonlinearprogram';



end