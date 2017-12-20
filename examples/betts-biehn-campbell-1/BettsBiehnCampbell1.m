%--------------------------------------------------------------------------
% BettsBiehnCampbell1.m
% J. T. Betts, N. Biehn, and S. L. Campbell, Convergence of Nonconvergent 
% IRK Discretizations of Optimal Control Problems with State Inequality 
% Constraints," SIAM Journal on Scientific Computing, vol. 23, no. 6, 
% pp. 1981-2007, 2002. doi: 10.1137/S1064827500383044
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function varargout = BettsBiehnCampbell1(varargin)

% default parameters
opts.plotflag = 1; % create the plots
opts.saveflag = 0;
opts.Defectmethod = 'TR';
opts.Quadmethod = 'CTR';
opts.NType = 'ED';
p.nt = 11; % number of nodes

% if input arguments are provided
% BettsBiehnCampbell1(p,p.nt,opts,opts.Quadmethod,opts.Defectmethod,opts.NType)
if nargin >= 1
    p = varargin{1};
end
if nargin >= 2
    p.nt = varargin{2};
end
if nargin >= 3
    opts = varargin{3};
end
if nargin >= 4
    opts.Quadmethod = varargin{4};
end
if nargin >= 5
    opts.Defectmethod = varargin{5};
end
if nargin >= 6
    opts.NType = varargin{6};
end
if nargin > 6
    warning('too many input arguments...');
end

% set current file name and path
[mpath,mname] = fileparts(mfilename('fullpath'));
opts.mpath = mpath;
opts.mname = mname;

%% setup
p.t0 = 34/15; p.tf =4;

% system dynamics
A = [0 1; 0 0]; 
B = [0;1]; 

% Lagrange term
L(1).left = 1; % control variables
L(1).right = 1; % control variables
L(1).matrix = 1e-3; % 1e-3*u^2
L(2).left = 2; % state variables
L(2).right = 2; % state variables
L(2).matrix = [1 0; 0 0]; % y1^2

% initial values
LB(1).right = 4; % initial states
LB(1).matrix = [302399/50625; 70304/3375];
UB(1).right = 4; % initial states
UB(1).matrix = [302399/50625; 70304/3375];

% simple bound path constraint
LB(2).right = 2; % states
LB(2).matrix = {@(t) 15 - (t-4).^4;-Inf};

% combine
setup.A = A; setup.B = B; setup.L = L; 
setup.LB = LB; setup.UB = UB; setup.p = p;

%% solve
[T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = BettsBiehnCampbell1_output(T,U,Y,P,F,p,opts);
if nargout == 1
	varargout{1} = O;
end

%% plot
BettsBiehnCampbell1_plot(T,U,Y,P,F,p,opts,sol)