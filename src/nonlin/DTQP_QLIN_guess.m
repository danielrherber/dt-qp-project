%--------------------------------------------------------------------------
% DTQP_QLIN_guess.m
% Construct and solve the quasilinearization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,opts] = DTQP_QLIN_guess(setup,opts,o,D2)

% extract
n = setup.n; nu = n.nu; ny = n.ny; np = n.np;
nt = opts.dt.nt;

T = linspace(setup.t0,setup.tf,nt)';

% (potentially) initial guess values for multipliers
if isfield(opts.method,'sqpflag') && opts.method.sqpflag
    opts = guessMultipliers(opts,n,T,D2);
end

% (potentially) used-defined initial guess
if isfield(setup,'p') && isfield(setup.p,'guess')

    % interpolate initial guess matrix
    X0 = interp1([setup.t0 setup.tf],setup.p.guess,T);

    % extract
    U = X0(:,1:nu);
    Y = X0(:,nu+1:nu+ny);
    P = X0(1,nu+ny+1:end);

    return
end

% initial guess values for controls, states, and parameters
T = linspace(setup.t0,setup.tf,nt)';
U = ones(nt,nu);
Y = ones(nt,ny);
P = ones(np,1);

% TODO: add more initial guess options

end

% construct a guess multiplier structure
function opts = guessMultipliers(opts,o,T,D2)

% number of optimization variables
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