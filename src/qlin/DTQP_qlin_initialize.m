%--------------------------------------------------------------------------
% DTQP_qlin_initialize.m
% Construct the quasilinearization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup,opts] = DTQP_qlin_initialize(setup,opts)

% initialize
qlinflag = false; % quasilinearization not needed
lqdoflag = true; % LQDO problem elements only

% check if the symbolic field is present
if ~isfield(setup,'symb')
    setup(1).symb = [];
    opts.qlin.qlinflag = qlinflag;
    opts.qlin.lqdoflag = lqdoflag;
    return % no symbolic operations needed
else
    symb = setup.symb;
end

% extract
o = symb.o; % NOTE: is this always present?

% set default field values if not present
if ~isfield(o,'ny')
    o.ny = 0;
end
if ~isfield(o,'nu')
    o.nu = 0;
end
if ~isfield(o,'np')
    o.np = 0;
end
if ~isfield(symb,'param')
    symb.param = [];
end

% assign
symb.o = o;

% quadracize the objective term
if isfield(symb,'Ob')
    % quadraticization of a nonlinear scalar equation
    form = 4;

    % quadracize
    L = DTQP_qlin_taylor(symb.Ob,form,o);

    % assign
    symb.L = L;

    % quasilinearization code will be needed
    qlinflag = true;

    % check if the outputs are all numeric
    lqdoflag = lqdoflag && CheckCellNumeric(L.H,L.G,L.C);
end

% linearize the dynamics
if isfield(symb,'D')
    % linearization of a nonlinear vector of equations
    form = 3;

    % linearize
    Linf = DTQP_qlin_taylor(symb.D,form,o);

    % assign
    symb.Linf = Linf;

    % quasilinearization code will be needed
    qlinflag = true;

    % check if the outputs are all numeric
    lqdoflag = lqdoflag && CheckCellNumeric(Linf.A,Linf.B,Linf.G,Linf.d);
end

% (NOT VERIFIED) linearize the control constraints
if isfield(symb,'C')
    % linearization of a nonlinear vector of equations
    form = 3;

    % linearize
    Linc = DTQP_qlin_taylor(symb.C,form,o);

    % assign
    symb.c = Linc;

    % quasilinearization is needed
    qlinflag = true;
end

% (NOT VERIFIED) linearize the state constraints
if isfield(symb,'Y')
    % linearization of a nonlinear vector of equations
    form = 3;

    % linearize
    Liny = DTQP_qlin_taylor(symb.Yi,form,o);

    % assign
    symb.y = Liny;

    % quasilinearization is needed
    qlinflag = true;
end

% assign
setup.symb = symb;
opts.qlin.qlinflag = qlinflag;
opts.qlin.lqdoflag = lqdoflag;

end

% check if the cell matrix can be converted to a numeric matrix
function flag = CheckCellNumeric(varargin)

% try to convert each input
try
    flag = true; % all inputs are numeric
    for k = 1:nargin
        A = cell2mat(varargin{k});
        if ~isnumeric(A)
            flag = false; % some inputs are not numeric
            break
        end
    end
catch
    flag = false; % some inputs are not numeric
end

end