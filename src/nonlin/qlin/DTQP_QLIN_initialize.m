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
if ~isfield(setup,'symb')
    setup(1).symb = [];
    opts.method.qlinflag = false;
    opts.method.lqdoflag = true;
    return % no symbolic operations needed
else
    symb = setup.symb;
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
if ~isfield(symb,'paramstr')
    n.paramstr = [];
    n.param = [];
else
    n.paramstr = symb.paramstr;
    n.param = symb.param;
end

% assign
setup.n = n;

% initialize the problem data for the selected method
switch opts.method.form
    %----------------------------------------------------------------------
    case 'qlin'
    % quasilinearization
    [setup.symb,opts] = DTQP_QLIN_initialize_qlin(symb,n,opts);
    %----------------------------------------------------------------------
    case 'fmincon'
    % fmincon nonlinear solver
%     [setup.symb,opts] = DTQP_QLIN_initialize_fmincon(symb,nopts);
end

% (potentially) end the timer
if (opts.general.displevel > 0) % minimal
    opts.timer.sym = opts.timer.sym + toc(opts.timer.t2); % start timer
end

end

function [symb,opts] = DTQP_QLIN_initialize_qlin(symb,n,opts)

% extract
sqpflag = opts.method.sqpflag;

% initialize
qlinflag = false; % quasilinearization not needed
lqdoflag = true; % LQDO problem elements only

% quadracize the objective term
if isfield(symb,'Ob')

    % quadraticization of a nonlinear scalar equation
    form = 4;

    % quadracize
    L = DTQP_QLIN_symb(symb.Ob,form,n,sqpflag);

    % assign
    symb.L = L;

    % quasilinearization code will be needed
    qlinflag = true;

    % check if the outputs are all numeric
    lqdoflag = lqdoflag && DTQP_checkCellNumeric(1,L.H,L.G,L.C);
end

% linearize the dynamics
if isfield(symb,'D')

    % linearization of a nonlinear vector of equations
    form = 3;

    % linearize
    Linf = DTQP_QLIN_symb(symb.D,form,n,sqpflag);

    % assign
    symb.Linf = Linf;

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