%--------------------------------------------------------------------------
% DTQP_solve.m
% Construct and solve a DO problem with DTQP
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

% (potentially) display banner
if (displevel > 1) % minimal

    disp('----------------------------------------------------------------')
    disp("                        █▀▄ ▀█▀ █▀█ █▀█                         ")
    disp("                        █▄▀  █  ▀▀█ █▀▀                         ")
    disp('Link: <a href = "https://github.com/danielrherber/dt-qp-project">https://github.com/danielrherber/dt-qp-project</a>')
    str = strings(0);
    str(1) = strcat("Machine: ",string(java.net.InetAddress.getLocalHost.getHostName));
    str(2) = strcat("Date: ",datestr(now,'yyyy/mm/dd HH:MM'));
    str = strjoin(str,' | ');
    disp(str)
    disp('----------------------------------------------------------------')

end

% (potentially) start the timers
if (displevel > 0) % minimal

    opts.timer.t1 = tic; % start timer
    opts.timer.sym = 0;
    opts.timer.create = 0;
    opts.timer.qpsolver = 0;
    opts.timer.t3 = tic; % start timer

end

% solve the problem
[T,U,Y,P,F,in,opts] = DTQP_MESH(setup,opts);

% (potentially) end the timers
if (displevel > 0) % minimal

    opts.timer.create = opts.timer.create + toc(opts.timer.t3); % add
    opts.timer.total = toc(opts.timer.t1); % start timer
    in(1).QPcreatetime = opts.timer.sym + opts.timer.create;
    in(1).QPsolvetime = opts.timer.qpsolver;

end

% (potentially) display to the command window
if (displevel > 1) % verbose

    disp(['█ Total time: ', num2str(opts.timer.total), ' s'])
    disp('----------------------------------------------------------------')

end

end