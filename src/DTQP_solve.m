%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve the LQDO problem
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,p,opts] = DTQP_solve(setup,opts)

    % check if it is a multiphase problem
    if length(setup) > 1
        % construct and solve the multiphase problem
        [T,U,Y,P,F,p,opts] = DTQP_multiphase(setup,opts);
    else % single phase

        % intialize some stuff
        [setup,opts] = DTQP_default_opts(setup,opts);
        
        % transcribe the problem
        [H,f,c,A,b,Aeq,beq,lb,ub,setup,p,opts] = DTQP_create(setup,opts);

        % (optional) simple scaling
        if ~isempty(setup.scaling)
            [H,f,A,b,Aeq,beq,lb,ub,p,s] = DTQP_scaling(H,f,A,b,Aeq,beq,lb,ub,p,setup.scaling);
        end

        % (optional) reordering of optimization variables 
        if opts.reorder
            [H,f,c,A,b,Aeq,beq,lb,ub,~,~] = DTQP_reorder(H,f,c,A,b,Aeq,beq,lb,ub,p,[],0);
        end

        % solve the optimization problem
        [X,F,opts] = DTQP_solver(H,f,A,b,Aeq,beq,lb,ub,opts);
        
        % (optional) restore ordering
        if opts.reorder
            [~,~,~,~,~,~,~,~,~,~,X] = DTQP_reorder([],[],[],[],[],[],[],[],[],p,X,1);
        end

        % optional unscale solution
        if ~isempty(setup.scaling)
            X = X.*s; % unscale optimization variables
        end    

        % add the constant term to objective function
        F = F + c;        

        % return optimal controls, states, and parameters
        T = p.t;
        U = reshape(X(1:p.nu*p.nt),p.nt,p.nu); % controls
        Y = reshape(X(p.nu*p.nt+1:(p.nu+p.ns)*p.nt),p.nt,p.ns); % states
        P = reshape(X((p.nu+p.ns)*p.nt+1:(p.nu+p.ns)*p.nt+p.np),p.np,1); % parameters
        
        % check for zero-order hold method and nan final controls
        if strcmp(opts.Defectmethod,'ZO') || strcmp(opts.Defectmethod,'EF')
            U(end,:) = nan;
        end
    end
    
end