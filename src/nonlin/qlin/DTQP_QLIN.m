%--------------------------------------------------------------------------
% DTQP_QLIN.m
% Construct and solve the quasilinearization problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_QLIN(setup,opts)

% initialize some stuff for quasilinearization
[setup,opts] = DTQP_QLIN_initialize(setup,opts);

% extract
displevel = opts.general.displevel;
plotflag = opts.general.plotflag;
lqdoflag = opts.method.lqdoflag;
tolerance = opts.method.tolerance;
improveX0flag = opts.method.improveguess;
deltascaleflag = opts.method.deltascaleflag;
sqpflag = opts.method.sqpflag;
imax = opts.method.maxiters;
symb = setup.symb;
o = setup.n;
param = o.param;

% check if this is an lqdo problem
if lqdoflag
    imax = 0; % only the initial iteration needed
end

% linearized dynamics
D2 = [];
if isfield(symb,'Linf')
    Dflag = true;
    D = symb.Linf;
    DA = D.A; DB = D.B; DG = D.G; Dd = D.d;
    if sqpflag
       D2 = D.D2;
    end
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
[T,U,Y,P,opts] = DTQP_QLIN_guess(setup,opts,o,D2);

% initialize
iter = 0;
Fold = 0;
F = inf;
opts.reduction = inf;

% potentially create time-varying parameter matrix
if isa(param,'cell')
    param = squeeze(DTQP_tmatrix(param,setup.p,T));
end

% quasilinearization
while (tolerance <= abs(F-Fold)) && (iter <= imax)

    % store previous F value
    Fold = F;

    % copy setup
    setupi = setup;

    if iter == 0
    	if improveX0flag
	        [U,Y,P,~] = DTQP_QLIN_improve_initial_point(setupi,opts,T,U,Y,P,param,Dflag,DA,DB,DG,Dd);
    	end
    end

    % construct previous solution vector
    Pe = repelem(P(:)',opts.dt.nt,1);
    X = [U,Y,Pe];

    % update dynamics based on previous solution vector
    if Dflag
        setupi = DTQP_QLIN_update_dynamics(setupi,DA,DB,DG,Dd,T,X,param);

        % update second derivative matrix for state derivative function
        if sqpflag
            % NEED
            D2s = cell(size(D2));
            for i = 1:length(D2)
                D2s{i} = DTQP_QLIN_update_tmatrix(D2{i},T,X,param);
            end
            setupi.D2 = D2s;

        end

    end

    % update Lagrange terms
    if Lflag
        setupi = DTQP_QLIN_update_lagrange(setupi,LH,LG,LC,o,T,X,param);
    end

    % TODO: update Mayer terms

    % TODO: update all other constraints
    % setup = DTQP_QLIN_update_control_constraint(setup,opts);
    % setup = DTQP_QLIN_update_state_constraint(setup,opts);

    % (potentially) shift optimization variables by previous solution
    if deltascaleflag
        setupi.scaling(1).right = 1; % controls
        setupi.scaling(1).constant = U;
        setupi.scaling(2).right = 2; % states
        setupi.scaling(2).constant = Y;
        setupi.scaling(3).right = 3; % parameters
        setupi.scaling(3).constant = P;
    end

    % solve the LQDO problem using DT and (potentially) mesh refinement
    [T,U,Y,P,F,in,opts] = DTQP_MESH(setupi,opts);

    % check if the previous problem failed
    if isnan(F)
        if (displevel > 0) % minimal
            disp("WARNING: did not solve with given opts.solver.tolerance")
        end

        % extract
        qptolerance = opts.solver.tolerance;

        % new "loose" tolerance
        opts.solver.tolerance = 1e-4;

        % solve the LQDO problem using DT and (potentially) mesh refinement
        [T,U,Y,P,F,in,opts] = DTQP_MESH(setupi,opts);

        % reassign
        opts.solver.tolerance = qptolerance;

    end

    % (potentially) plot current iteration
    if (plotflag > 0)
        DTQP_QLIN_plots(T,Y,U,P,iter+1)
    end

    % (potentially) display to  command window
    if (displevel > 1) % minimal
        qlinDispFun(iter,F,abs(F-Fold),in)
    end

    % increment iteration counter
    iter = iter + 1;

    opts.reduction = abs(F-Fold);
end

% add outputs
in.output.iterations = iter;
in.output.funcCount = nan;

% (potentially) plot final iteration
if (plotflag > 0)
    DTQP_QLIN_plots(T,Y,U,P,-iter)
end

end

function qlinDispFun(iter,F,E,in)

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
if isfield(in.output,'sqppenalty')
    Pstr = sprintf('%1.3e',in.output.sqppenalty);
else
    Pstr = blanks(length(sprintf('%1.3e',1)));
end

% initial headers
if iter == 0
    disp('-------------- qlin progress ---------------')
    disp('| iter |         F |        dF |  sqp pen  |')
end

% values
disp(['| ',sprintf(' %3i',iter),...
    ' | ', Fstr,...
    ' | ', Estr,...
    ' | ', Pstr,...
    ' | ']);

end