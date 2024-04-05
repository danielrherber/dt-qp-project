%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve a dynamic optimization (DO) problem with DTQP
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts)

% (potentially) convert input structure
setup = DTQP_setup_convert(setup);

% initialize options and other things
[setup,opts] = DTQP_default_opts(setup,opts);

% extract
displevel = opts.general.displevel;

% (potentially) display banner
if (displevel > 1) % verbose
    flag = 'line'; DTQP_commandWindowTasks %#ok<NASGU>
    flag = 'banner'; DTQP_commandWindowTasks %#ok<NASGU>
    flag = 'link'; DTQP_commandWindowTasks %#ok<NASGU>
    flag = 'info'; DTQP_commandWindowTasks %#ok<NASGU>
    flag = 'line'; DTQP_commandWindowTasks %#ok<NASGU>
end

% (potentially) start the timers
if (displevel > 0) % minimal
    opts.timer.t1 = tic; % start timer (for total time)
    opts.timer.sym = 0;
    opts.timer.create = 0;
    opts.timer.qpsolver = 0;
    opts.timer.t3 = tic; % start timer (for intermediate times)
end

% solve the problem
[T,U,Y,P,F,in,opts] = DTQP_MESH(setup,opts);

% (potentially) end the timers
if (displevel > 0) % minimal
    opts.timer.create = opts.timer.create + toc(opts.timer.t3); % add
    opts.timer.total = toc(opts.timer.t1); % end timer
    in.QPcreatetime = opts.timer.sym + opts.timer.create; % add
    in.QPsolvetime = opts.timer.qpsolver;
end

% (potentially) display to the command window
if (displevel > 1) % verbose
    flag = 'total-time'; DTQP_commandWindowTasks %#ok<NASGU>
    flag = 'line'; DTQP_commandWindowTasks %#ok<NASGU>
end

end