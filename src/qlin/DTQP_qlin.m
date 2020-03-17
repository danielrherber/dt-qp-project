%--------------------------------------------------------------------------
% DTQP_qlin.m
% Construct and solve the quasilinearization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_qlin(setup,opts)

% get default options
opts = DTQP_qlin_default_opts(opts);

% extract
displevel = opts.general.displevel;
plotflag = opts.general.plotflag;
lqdoflag = opts.qlin.lqdoflag;
tolerance = opts.qlin.tolerance;
imax = opts.qlin.imax;
symb = setup.symb;
param = symb.param;
o = symb.o;

% check if this is an lqdo problem
if lqdoflag
    imax = 0; % only the initial iteration needed
end

% linearized dynamics
if isfield(symb,'Linf')
    Dflag = true;
    D = symb.Linf;
    DA = D.A; DB = D.B; DG = D.G; Dd = D.d;
else
    Dflag = false;
end

% quadracized Lagrange term
if isfield(symb,'L')
    Lflag = true;
    L = symb.L;
    LH = L.H; LG = L.G; LC = L.C;
else
    Lflag = false;
end

% TODO: extract other constraints

% initial guess values for controls, states, and parameters
[T,U,Y,P] = DTQP_qlin_guess(setup,opts,o);

% initialize
iter = 0;
Fold = 0;
F = inf;

% quasilinearization
while (tolerance <= abs(F-Fold)) && (iter <= imax)

    % store previous F value
    Fold = F;

    % copy setup
    setupi = setup;

    % construct previous solution vector
    P = repelem(P',opts.dt.nt,1);
    X = [U,Y,P];

    % update dynamics based on previous solution vector
    if Dflag
        setupi = DTQP_qlin_updateDynamics(setupi,DA,DB,DG,Dd,T,X,param);
    end

    % update Lagrange terms
    if Lflag
        setupi = DTQP_qlin_updateLagrange(setupi,LH,LG,LC,o,T,X,param);
    end

    % TODO: update Mayer terms

    % TODO: update all other constraints
    % setup = DTQP_qlin_updateControlConstraint(setup,opts);
    % setup = DTQP_qlin_updateStateConstraint(setup,opts);

    % solve the LQDO problem using DT and (potentially) mesh refinement
    [T,U,Y,P,F,in,opts] = DTQP_meshr(setupi,opts);

    % (potentially) plot current iteration
    if (plotflag > 0)
        DTQP_qlin_plots(T,Y,U,P,iter+1)
    end

    % (potentially) display to  command window
    if (displevel > 0) % minimal
        qlinDispFun(iter,F,abs(F-Fold))
    end

    % increment iteration counter
    iter = iter + 1;

end

% (potentially) plot final iteration
if (plotflag > 0)
    DTQP_qlin_plots(T,Y,U,P,0)
end

end

function qlinDispFun(iter,F,E)

% handle edge cases
if isinf(E)
    Estr = blanks(length(sprintf('%1.3e',1)));
else
    Estr = sprintf('%1.3e',E);
end
if isinf(F)
    Fstr = blanks(length(sprintf('%1.3e',1)));
else
    Fstr = sprintf('%1.3e',F);
end

% initial headers
if iter == 0
    disp('-------- qlin progress ---------')
    disp('| iter |         F |        dF |')
end

% values
disp(['| ',sprintf(' %3i',iter),...
    ' | ', Fstr,...
    ' | ', Estr,...
    ' | ']);

end