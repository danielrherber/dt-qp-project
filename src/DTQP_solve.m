%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve a (LQDO) problem with DTQP
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts)

% initialize some stuff
[setup,opts] = DTQP_default_opts(setup,opts);

% extract
displevel = opts.general.displevel;

% (potentially) start the timers
if (displevel > 0) % minimal
    opts.timer.t1 = tic; % start timer
    opts.timer.sym = 0;
    opts.timer.create = 0;
    opts.timer.qpsolver = 0;
end

% initialize some stuff for quasilinearization
[setup,opts] = DTQP_qlin_initialize(setup,opts);

% (potentially) start the timer
if (opts.general.displevel > 0) % minimal
    opts.timer.t3 = tic; % start timer
end

% check if quasilinearization is needed
if opts.qlin.qlinflag
    % solve the problem with quasilinearization
    [T,U,Y,P,F,in,opts] = DTQP_qlin(setup,opts);
else
    % solve the LQDO problem (potentially) using mesh refinement
    [T,U,Y,P,F,in,opts] = DTQP_meshr(setup,opts);
end

% (potentially) end the timers
if (displevel > 0) % minimal
    opts.timer.create = opts.timer.create + toc(opts.timer.t3); % add
    opts.timer.total = toc(opts.timer.t1); % start timer
    opts.timer = rmfield(opts.timer,{'t1','t2','t3'}); % remove fields
    in(1).QPcreatetime = opts.timer.sym + opts.timer.create;
    in(1).QPsolvetime = opts.timer.qpsolver;
end

% (potentially) display to the command window
if (displevel > 1) % verbose
    disp('----------------------------------------')
    disp(['total time: ', num2str(opts.timer.total), ' s'])
end

% disp(opts.timer.total)
% disp(opts.timer.sym + opts.timer.create + opts.timer.qpsolver)

end