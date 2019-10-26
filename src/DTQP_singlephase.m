%--------------------------------------------------------------------------
% DTQP_singlephase.m
% Construct and solve the LQDO problem for a given mesh (single phase)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_singlephase(setup,opts)

    % transcribe the problem
    [H,f,c,A,b,Aeq,beq,lb,ub,setup,in,opts] = DTQP_create(setup,opts);

    % (optional) simple scaling
    if ~isempty(setup.scaling)
        [H,f,A,b,Aeq,beq,lb,ub,in,s] = DTQP_scaling(H,f,A,b,Aeq,beq,lb,ub,in,setup.scaling);
    end

    % (optional) reordering of optimization variables 
    if opts.qp.reorder
        [H,f,c,A,b,Aeq,beq,lb,ub] = DTQP_reorder(in,H,f,c,A,b,Aeq,beq,lb,ub);
    end

    % solve the optimization problem
    [X,F,in,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,in,opts);

    % (optional) restore ordering
    if opts.qp.reorder
        X = DTQP_reorder(in,X);
    end

    % optional unscale solution
    if ~isempty(setup.scaling)
        X = X.*s; % unscale optimization variables
    end    

    % add the constant term to objective function
    F = F + c;        

    % return optimal controls, states, and parameters
    T = in.t;
    U = reshape(X(1:in.nu*in.nt),in.nt,in.nu); % controls
    Y = reshape(X(in.nu*in.nt+1:(in.nu+in.ny)*in.nt),in.nt,in.ny); % states
    P = reshape(X((in.nu+in.ny)*in.nt+1:(in.nu+in.ny)*in.nt+in.np),in.np,1); % parameters

    % check for zero-order hold method and nan final controls
    if strcmpi(opts.dt.defects,'ZO') || strcmpi(opts.dt.defects,'EF')
        U(end,:) = nan;
    end

end