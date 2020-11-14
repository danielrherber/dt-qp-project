%--------------------------------------------------------------------------
% DTQP_TEST_bench_LQDO.m
% Benchmark code using the LQDO problems
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clear; clc; close all;

optsnum = 1; % options test (see below)
testnum = 3; % test scenario (see below)

% options setup
switch optsnum
    %------------------------------------------------------------------
    case 1
    opts.dt.defects = 'TR';
    opts.dt.quadrature = 'CTR';
    opts.dt.mesh = 'ED';
    %------------------------------------------------------------------
    case 2
    opts.dt.defects = 'HS';
    opts.dt.quadrature = 'CQHS';
    opts.dt.mesh = 'ED';
    %------------------------------------------------------------------
    case 3
    opts.dt.defects = 'ZO';
    opts.dt.quadrature = 'CEF';
    opts.dt.mesh = 'ED';
    %------------------------------------------------------------------
    case 4
    opts.dt.defects = 'PS';
    opts.dt.quadrature = 'G';
    opts.dt.mesh = 'LGL';
end

% select solver
opts.solver.function = 'QUADPROG';
% opts.solver.function = 'QPOASES';
% opts.solver.function = 'CVX';

% set some other options
opts.general.plotflag = 0; % no plots
opts.general.displevel = 1; % only calculate compute times

% initialize
p = [];

% create the test scenario
switch testnum
    %----------------------------------------------------------------------
    case 1 % single function, multiple n values
    F = strings(1,0); N = cell(1,0);
    F(end+1) = "BrysonDenham";
    N{end+1} = [20,50,100];
    %----------------------------------------------------------------------
    case 2 % several functions, single n values
    F = strings(1,0); N = cell(1,0);
    F(end+1) = "AndersonMoore64";
    N{end+1} = [100,1000,10000];
    F(end+1) = "BettsBiehnCampbell1";
    N{end+1} = [100,1000,10000];
    %----------------------------------------------------------------------
    case 3 % primary test (only some problems in the benchmark suite)
    F = strings(1,0); N = cell(1,0);
    F(end+1) = "BettsBiehnCampbell1";
    N{end+1} = [100,1000,10000];
    F(end+1) = "BrysonDenham";
    N{end+1} = [100,1000,10000];
    F(end+1) = "BrysonHo116";
    N{end+1} = [100,1000,10000];
    F(end+1) = "DTQP1";
    N{end+1} = [100,1000,10000];
    F(end+1) = "GasAbsorber";
    N{end+1} = [100,1000,10000];
    F(end+1) = "HDAE";
    N{end+1} = [100,1000,10000];
    F(end+1) = "LQRInhomogeneous";
    N{end+1} = [100,1000,10000];
    F(end+1) = "LQRstandard";
    N{end+1} = [100,1000,10000];
    F(end+1) = "MinimumEnergyTransfer";
    N{end+1} = [100,1000,10000];
    F(end+1) = "MultiphaseParameter";
    N{end+1} = [100,1000,10000];
    F(end+1) = "Nagurka";
    N{end+1} = [100,1000,10000];
    F(end+1) = "OutputTracking";
    N{end+1} = [100,1000,10000];
    %----------------------------------------------------------------------
end

% n format string
idxFormat = ['%0',num2str(max(1,ceil(log10(max(cellfun(@max, N)))))),'i']; % pad with correct number of zeros

% initialize incorrect counter
fcounter = 0;

% go through each function to test
for idx = 1:length(F)

    % go through each number of time points
    for n = 1:length(N{idx})

        % try to run the test
        try

            % number of time points
            opts.dt.nt = N{idx}(n);

            % reduce number of nodes with the PS method (heuristic formula)
            if strcmpi(opts.dt.defects,"PS")
                opts.dt.nt = round(1.034*opts.dt.nt^0.574 - 4.545);
            end

            % run test
            eval(strcat("o = ",F(idx),"(p,opts);"));

            % find times
            QPcreatetime = o(strcmpi({o.label},'QPcreatetime')).value;
            QPsolvetime = o(strcmpi({o.label},'QPsolvetime')).value;

            % error in the objective (if available)
            IF = strcmpi({o.label},'F');

            % check if a solution was found
            if any(IF)
                Ferror = o(IF).value;
            else
                Ferror = "-";
            end

            % test passed
            f = true;

        catch
            % failed
            f = false; QPcreatetime = -1; QPsolvetime = -1; F = " ";
        end

        % construct string
        str = strings(0);
        str(1) = F(idx); % sequence
        str(2) = string(opts.dt.nt);
        str(3) = num2str(QPcreatetime, '%10.3e'); % time
        str(4) = num2str(QPsolvetime, '%10.3e'); % time
        str(5) = num2str(Ferror, '%10.2e'); % F error
        str(6) = strcat(opts.dt.defects,"-",opts.dt.quadrature,"-",opts.dt.mesh);
        str(7) = string(java.net.InetAddress.getLocalHost.getHostName);
        str(8) = datestr(now,'yyyy/mm/dd');
        str = strjoin(str,' : ');

        % only if an error
        if ~f
            str = strcat(str," <- NOT CORRECT!!!" );
            fcounter = fcounter + 1;
        end

        % print
        disp(str)

    end
end

% print
disp(strcat(string(fcounter)," tests reported incorrect values"))