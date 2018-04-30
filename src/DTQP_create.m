%--------------------------------------------------------------------------
% DTQP_create.m
% Create the matrices that represent the quadratic program (QP)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,c,A,b,Aeq,beq,lb,ub,setup,p,opts] = DTQP_create(setup,opts)
    
    % extract parameter structure
    p = setup.p;

    % potentially start the timer
    if (opts.general.displevel > 0) % minimal
        tic % start timer
    end

    % initialize some stuff
    [setup,p] = DTQP_initialize(setup,opts.dt);
    
    % objective function (maximum quadratic terms)
    H = DTQP_createH(setup.L,setup.M,p,opts); % create Hessian
    f = DTQP_createf(setup.l,setup.m,p,opts); % create gradient
    c = DTQP_createc(setup.cL,setup.cM,p,opts); % determine constants

    % constraints (maximum linear terms)   
    % create defect constraints
    [Aeq1,beq1] = DTQP_defects(setup.A,setup.B,setup.G,setup.d,p,opts);
	% create linear path and boundary equality constraints
    [Aeq2,beq2] = DTQP_create_YZ(setup.Y,p);
    % combine linear equality constraints
    Aeq = [Aeq1;Aeq2];
    beq = [beq1;beq2]; % Aeq*X = beq
    
	% create linear path and boundary inequality constraints
    [A,b] = DTQP_create_YZ(setup.Z,p); % A*X <= b
    
    % create simple bounds (box bounds)
    [lb,ub] = DTQP_create_bnds(setup.LB,setup.UB,p);
    
    % end the timer
    if (opts.general.displevel > 0) % minimal
        opts.QPcreatetime = toc;
    end
    
    % display to the command window
    if (opts.general.displevel > 1) % verbose
        disp(['QP creation time: ', num2str(opts.QPcreatetime), ' s'])
    end

end