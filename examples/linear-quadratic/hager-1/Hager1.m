%--------------------------------------------------------------------------
% Hager1.m
% W. W. Hager, Runge-Kutta Methods in Optimal Control and the Transformed
% Adjoint System," Numerische Mathematik, vol. 87, no. 2, pp. 247-282, 
% Dec. 2000. doi: 10.1007/s002110000178
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = Hager1(varargin)
% input arguments can be provided in the format 'Hager1(p,opts)'

% set local functions
ex_opts = @Hager1_opts; % options function
ex_output = @Hager1_output; % output function
ex_plot = @Hager1_plot; % plot function

% set p and opts (see local_opts)
[p,opts] = DTQP_standardizedinputs(ex_opts,varargin);

%% setup
t0 = 0; tf = 1; % time horizon

% system dynamics
A = 0.5; 
B = 1; 

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1; 
L(2).left = 2; % state variables
L(2).right = 1; % control variables
L(2).matrix = 1; 
L(3).left = 2; % state variables
L(3).right = 2; % state variables
L(3).matrix = 5/4; 

% initial conditions
LB(1).right = 4; % initial states
LB(1).matrix = 1;
UB(1).right = 4; % initial states
UB(1).matrix = 1;

% combine
setup.A = A; setup.B = B; setup.L = L;
setup.LB = LB; setup.UB = UB; setup.t0 = t0; setup.tf = tf; setup.p = p;

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
function opts = Hager1_opts
% test number
num = 1;

switch num
case 1
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 10;
case 2
    opts.dt.quadrature = 'CQHS';
    opts.dt.defects = 'HS';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
case 3
    opts.dt.quadrature = 'CTR';
    opts.dt.defects = 'TR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10;
end

end