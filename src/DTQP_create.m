%--------------------------------------------------------------------------
% DTQP_create.m
% Create the matrices that represent the quadratic program (QP)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [H,f,c,A,b,Aeq,beq,lb,ub,setup,in,opts] = DTQP_create(setup,opts)

    % potentially start the timer
    if (opts.general.displevel > 0) % minimal
        tstart = toc; % start timer
    end

    % initialize some stuff
    [setup,in] = DTQP_initialize(setup,opts.dt);
    
    % objective function (maximum quadratic terms)
    H = DTQP_createH(setup.L,setup.M,in,opts); % create Hessian
    f = DTQP_createf(setup.l,setup.m,in,opts); % create gradient
    c = DTQP_createc(setup.cL,setup.cM,in,opts); % determine constants

    % constraints (maximum linear terms)   
    % create defect constraints
    [Aeq1,beq1] = DTQP_defects(setup.A,setup.B,setup.G,setup.d,in,opts);
	% create linear path and boundary equality constraints
    [Aeq2,beq2] = DTQP_create_YZ(setup.Y,in);
    % combine linear equality constraints
    Aeq = [Aeq1;Aeq2];
    beq = [beq1;beq2]; % Aeq*X = beq
    
	% create linear path and boundary inequality constraints
    [A,b] = DTQP_create_YZ(setup.Z,in); % A*X <= b
    
    % create simple bounds (box bounds)
    [lb,ub] = DTQP_create_bnds(setup.LB,setup.UB,in);
    
    % end the timer
    if (opts.general.displevel > 0) % minimal
        in.QPcreatetime = toc - tstart;
    end
    
    % display to the command window
    if (opts.general.displevel > 1) % verbose
        disp(['QP creation time: ', num2str(in.QPcreatetime), ' s'])
    end

end