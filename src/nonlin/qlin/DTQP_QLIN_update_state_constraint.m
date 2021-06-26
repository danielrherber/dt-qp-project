%--------------------------------------------------------------------------
% DTQP_QLIN_update_state_constraint.m
% Convert the nonlinear state constraint to a Linear constraint for the
% LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup] = DTQP_QLIN_update_state_constraint(setup,opts)

% get the respective Nonlinear variables
element = setup.element;
Liny = element.y;
oy = element.oy;

% get the linearized matrix values
B = rand(opts.dt.nt,oy.ny);
T = linspace(setup.auxdata.t0,setup.auxdata.tf,opts.dt.nt)';
X = [B];
param = [];
A = TaylorConvert(Liny.A,T,X,param);
d = TaylorConvert(Liny.d,T,X,param);

% Linear equality constraints
Y(2).linear(1).right = 5; Y(2).linear(1).matrix = A;
Y(2).b = d;

% return the values
setup.Y = Y;

end