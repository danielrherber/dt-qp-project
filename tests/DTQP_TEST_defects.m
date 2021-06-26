%--------------------------------------------------------------------------
% DTQP_TEST_defects.m
% Test defect constraint creation and error
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:7; % all
% tests = 1:5; % single-step methods
% tests = 6:7; % pseudospectral methods

% number of time points
opts.dt.nt = 10;

% problem structure
[setup,opts,F,auxdata,T] = problem(opts);

% initialize
l1str{1} = '$\xi_{actual}$';
l2str = cell(0);
c = lines(length(tests));

% go through the tests
for k = 1:length(tests)

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        opts.dt.defects = 'ZO';
        opts.dt.mesh = 'ED';
        l1str{end+1} = '$\xi_{ZO}$';
        l2str{end+1} = '$E_{ZO}$';
        %------------------------------------------------------------------
        case 2
        opts.dt.defects = 'EF';
        opts.dt.mesh = 'ED';
        l1str{end+1} = '$\xi_{EF}$';
        l2str{end+1} = '$E_{EF}$';
        %------------------------------------------------------------------
        case 3
        opts.dt.defects = 'TR';
        opts.dt.mesh = 'ED';
        l1str{end+1} = '$\xi_{TR}$';
        l2str{end+1} = '$E_{TR}$';
        %------------------------------------------------------------------
        case 4
        opts.dt.defects = 'HS';
        opts.dt.mesh = 'ED';
        l1str{end+1} = '$\xi_{HS}$';
        l2str{end+1} = '$E_{HS}$';
        %------------------------------------------------------------------
        case 5
        opts.dt.defects = 'RK4';
        opts.dt.mesh = 'ED';
        l1str{end+1} = '$\xi_{RK4}$';
        l2str{end+1} = '$E_{RK4}$';
        %------------------------------------------------------------------
        case 6
        opts.dt.defects = 'PS';
        opts.dt.mesh = 'LGL';
        l1str{end+1} = '$\xi_{PS-LGL}$';
        l2str{end+1} = '$E_{PS-LGL}$';
        %------------------------------------------------------------------
        case 7
        opts.dt.defects = 'PS';
        opts.dt.mesh = 'CGL';
        l1str{end+1} = '$\xi_{PS-CGL}$';
        l2str{end+1} = '$E_{PS-CGL}$';
    end

    % run the test and time
    [t{k},~,X{k},~,~,~,~] = DTQP_solve(setup,opts);

    % test analysis
    figure(1)
    plot(t{k},X{k},'linewidth',1.5,'Color',c(k,:));

    figure(2); hold on
    semilogy(t{k},abs(X{k} - F(t{k},auxdata.x0))+1e-20,...
        'linewidth',1.5,'Color',c(k,:));
end

figure(1)
hl = legend(l1str,'location','Best');
hl.FontSize = 12; % change legend font size
hl.EdgeColor = 'k'; % change the legend border to black (not a dark grey)

figure(2)
hl = legend(l2str,'location','Best');
hl.FontSize = 12; % change legend font size
hl.EdgeColor = 'k'; % change the legend border to black (not a dark grey)

function [setup,opts,F,auxdata,T] = problem(opts)

% symbolically generate the correct function
syms t x0
syms x(t)
f = t*sin(t);
F = dsolve( diff(x,t) == f, x(0) == x0);
f = matlabFunction(f);
F = matlabFunction(F);

% tunable parameters
setup.t0 = 0; setup.tf = 10; % time horizon
auxdata.x0 = 1;
setup.auxdata = auxdata;

% system dynamics
setup.d{1,1} = f;

% initial state
setup.Y(1).linear.right = 4; % initial states
setup.Y(1).linear.matrix = 1;
setup.Y(1).b = auxdata.x0;

% dummy objective
setup.M(1).left = 5;
setup.M(1).right = 5;
setup.M(1).matrix = 1;

% don't display anything
opts.general.displevel = 0;

bcolor = [0 0 0]; % black color
fontlabel = 20; % x,y label font size
fonttick = 12; % x,y tick font size

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

%--- state figure
hf = figure(1); hold on % create a new figure and save handle
hf.Color = [1 1 1]; % change the figure background color
hf.Position = [200 200 550 400]; % set figure size and position

% high resolution time vector
T = linspace(setup.t0,setup.tf,10000);

% plot actual solution
plot(T,F(T,auxdata.x0),'k','linewidth',2); hold on

% post tasks
xlabel('$t$') % create x label
ylabel('$\xi$') % create y label
ha = gca; % get current axis handle
ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
ha.YAxis.Color = bcolor; % change the y axis color to black (not a dark grey)
ha.XAxis.FontSize = fonttick; % change x tick font size
ha.YAxis.FontSize = fonttick; % change y tick font size
ha.XAxis.Label.FontSize = fontlabel; % change x label font size
ha.YAxis.Label.FontSize = fontlabel; % change y label font size
ha.Layer = 'top'; % place the axes on top of the data

%--- error figure
hf = figure(2); hold on % create a new figure and save handle
hf.Color = [1 1 1]; % change the figure background color
hf.Position = [200 200 550 400]; % set figure size and position

% post tasks
xlabel('$t$') % create x label
ylabel('error') % create y label
ha = gca; % get current axis handle
ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
ha.YAxis.Color = bcolor; % change the y axis color to black (not a dark grey)
ha.XAxis.FontSize = fonttick; % change x tick font size
ha.YAxis.FontSize = fonttick; % change y tick font size
ha.XAxis.Label.FontSize = fontlabel; % change x label font size
ha.YAxis.Label.FontSize = fontlabel; % change y label font size
ha.Layer = 'top'; % place the axes on top of the data
ha.YScale = 'log';
end