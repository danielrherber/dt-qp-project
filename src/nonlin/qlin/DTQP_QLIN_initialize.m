%--------------------------------------------------------------------------
% DTQP_QLIN_initialize.m
% Construct the quasilinearization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup,opts] = DTQP_QLIN_initialize(setup,opts)

% (potentially) start the timer
if (opts.general.displevel > 0) % minimal
    opts.timer.t2 = tic; % start timer
end

% check if the symbolic field is present
if ~isfield(setup,'element')
    setup(1).element = [];
    opts.method.qlinflag = false;
    opts.method.lqdoflag = true;
    return % no symbolic operations needed
else
    element = setup.element;
end

% set default field values if not present
if ~isfield(setup,'n')
    n = [];
else
    n = setup.n;
end
if ~isfield(n,'ny')
    n.ny = 0;
end
if ~isfield(n,'nu')
    n.nu = 0;
end
if ~isfield(n,'np')
    n.np = 0;
end
if ~isfield(element,'parameter_list')
    n.parameter_list = [];
    n.param = [];
else
    n.parameter_list = element.parameter_list;
    n.param = element.parameter_values;
end

% assign
setup.n = n;

% initialize the problem data for the selected method
switch opts.method.form
    %----------------------------------------------------------------------
    case 'qlin'
    % quasilinearization
    [setup.element,opts] = DTQP_QLIN_initialize_qlin(element,n,opts);
    %----------------------------------------------------------------------
    case 'fmincon'
    % fmincon nonlinear solver
%     [setup.element,opts] = DTQP_QLIN_initialize_fmincon(element,nopts);
end

% (potentially) end the timer
if (opts.general.displevel > 0) % minimal
    opts.timer.sym = opts.timer.sym + toc(opts.timer.t2); % start timer
end

end

function [element,opts] = DTQP_QLIN_initialize_qlin(element,n,opts)

% extract
sqpflag = opts.method.sqpflag;

% initialize
qlinflag = false; % quasilinearization not needed
lqdoflag = true; % LQDO problem elements only

% quadracize the objective term
if isfield(element,'lagrange') && ~isempty(element.lagrange)

    % quadraticization of a nonlinear scalar equation
    form = 4;

    % quadracize
    L = DTQP_QLIN_symb(element.lagrange,form,n,sqpflag);

    % assign
    element.L = L;

    % quasilinearization code will be needed
    qlinflag = true;

    % check if the outputs are all numeric
    lqdoflag = lqdoflag && DTQP_checkCellNumeric(1,L.H,L.G,L.C);
end

% linearize the dynamics
if isfield(element,'dynamics') && ~isempty(element.dynamics)

    % linearization of a nonlinear vector of equations
    form = 3;

    % linearize
    Linf = DTQP_QLIN_symb(element.dynamics,form,n,sqpflag);

    % assign
    element.Linf = Linf;

    % quasilinearization code will be needed
    qlinflag = true;

    % check if the outputs are all numeric
    lqdoflag = lqdoflag && DTQP_checkCellNumeric(1,Linf.A,Linf.B,Linf.G,Linf.d);
end

% TODO: linearize nonlinear equality constraints

% TODO: linearize nonlinear inequality constraints

% assign
opts.method.qlinflag = qlinflag;
opts.method.lqdoflag = lqdoflag;

end