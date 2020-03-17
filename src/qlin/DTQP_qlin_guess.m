%--------------------------------------------------------------------------
% DTQP_qlin_guess.m
% Construct and solve the quasilinearization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P] = DTQP_qlin_guess(setup,opts,o)

% initial guess values for controls, states, and parameters
T = linspace(setup.p.t0,setup.p.tf,opts.dt.nt)';
U = ones(opts.dt.nt,o.nu);
Y = ones(opts.dt.nt,o.ny);
P = ones(o.np,1);

% TODO: add more initial guess options

end