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
function [T,U,Y,P,opts] = DTQP_qlin_guess(setup,opts,o,D2)

% initial guess values for controls, states, and parameters
T = linspace(setup.t0,setup.tf,opts.dt.nt)';
U = ones(opts.dt.nt,o.nu);
Y = ones(opts.dt.nt,o.ny);
P = ones(o.np,1);

% (potentially) initial guess values for multipliers
if opts.qlin.sqpflag
    opts = guessMultipliers(opts,o,T,D2);
end

% TODO: add more initial guess options

end

% construct a guess multiplier structure
function opts = guessMultipliers(opts,o,T,D2)

% number of optimization variablesa
nx = opts.dt.nt*(o.nu + o.ny);

% determine the number of defect constraints
switch opts.dt.defects
    case 'PS'
        ndefect = length(T);
    otherwise
        ndefect = length(T)-1;
end

% equality constraints
eqlin1 = zeros(ndefect*length(D2),1); % defect constraints
eqlin2 = zeros(0,1); % need to update, general equality constraints

% combine
lambda.ineqlin = zeros(0,1); % need to update, general inequality constraints
lambda.eqlin = [eqlin1;eqlin2]; % equality constraints
lambda.lower = zeros(nx,1); % need to update, simple lower bounds
lambda.upper = zeros(nx,1); % need to update, simple upper bounds

% assign
opts.lambda = lambda;

end