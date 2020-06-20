%--------------------------------------------------------------------------
% DTQP_QLIN_improve_initial_point.m
% Improve the initial guess for the controls, states, and parameters using
% a feasibility problem with sequential linearization
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [U,Y,P,lambda] = DTQP_QLIN_improve_initial_point(setup,opts,T,U,Y,P,param,Dflag,DA,DB,DG,Dd)

% remove objective function terms and scaling
setup.M = []; setup.L = []; setup.D2 = []; setup.scaling = [];

% construct previous solution vector
Pe = repelem(P',opts.dt.nt,1);
X = [U,Y,Pe];

% deactivate qlin flags
trustregionflag = opts.qlin.trustregionflag;
sqpflag = opts.qlin.sqpflag;
deltascaleflag = opts.qlin.deltascaleflag;
opts.qlin.trustregionflag = false;
opts.qlin.sqpflag = false;
opts.qlin.deltascaleflag = false;

% initialize
iter = 0;
Xold = inf;

% stopping parameters
imax = 10; % maximum number of iterations
tolerance = 1e-10; % minimum change from previous iteration

while (tolerance <= norm(X-Xold,inf)) && (iter <= imax)
    % store previous solution
    Xold = X;

    % update dynamics based on previous solution vector
    if Dflag
        setup = DTQP_QLIN_update_dynamics(setup,DA,DB,DG,Dd,T,X,param);
    end

    % solve the feasibility problem
    [T,U,Y,P,F,in,opts] = DTQP_MESH(setup,opts);

    % terminate if the problem was not solved
    if isnan(F)
        break
    end

    % construct previous solution vector
    Pe = repelem(P',opts.dt.nt,1);
    X = [U,Y,Pe];

    % increment iteration counter
    iter = iter + 1;

    % for debugging
    % disp(norm(X-Xold,inf))
    % plot(T,Y,'k'); hold on
end

% plot(T,Y,'r','linewidth',2); hold on

% reassign
opts.qlin.trustregionflag = trustregionflag;
opts.qlin.sqpflag = sqpflag;
opts.qlin.deltascaleflag = deltascaleflag;

% return multipliers
lambda = opts.lambda;

end