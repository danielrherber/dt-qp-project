%--------------------------------------------------------------------------
% DTQPtest_defects.m
% Test defect constraint creation and error
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all
clear
clc

% symbolically generate the correct function
syms t x0
syms x(t)
f = t*sin(t);
F = dsolve( diff(x,t) == f, x(0) == x0);
f = matlabFunction(f);
F = matlabFunction(F);

% number of nodes
opts.dt.nt = 25;

% tunable parameters
p.t0 = 0; p.tf = 10; % time horizon
p.x0 = 1;
setup.p = p;

% system dynamics
setup.d{1,1} = f;

% initial state
setup.Y(1).linear.right = 4; % initial states
setup.Y(1).linear.matrix = 1;
setup.Y(1).b = p.x0;

% dummy objective
setup.M(1).left = 5;
setup.M(1).right = 5;
setup.M(1).matrix = 1;

% don't display anything
opts.general.displevel = 0;

%% ZO
opts.dt.defects = 'ZO';
opts.dt.mesh = 'ED';
[tzo,~,Xzo,~,~,~,~] = DTQP_solve(setup,opts);

%% EF
opts.dt.defects = 'EF';
opts.dt.mesh = 'ED';
[tef,~,Xef,~,~,~,~] = DTQP_solve(setup,opts);

%% TR
opts.dt.defects = 'TR';
opts.dt.mesh = 'ED';
[ttr,~,Xtr,~,~,~,~] = DTQP_solve(setup,opts);

%% HS
opts.dt.defects = 'HS';
opts.dt.mesh = 'ED';
[ths,~,Xhs,~,~,~,~] = DTQP_solve(setup,opts);

%% RK4
opts.dt.defects = 'RK4';
opts.dt.mesh = 'ED';
[trk4,~,Xrk4,~,~,~,~] = DTQP_solve(setup,opts);

%% PS-LGL
opts.dt.defects = 'PS';
opts.dt.mesh = 'LGL';
[tlgl,~,Xlgl,~,~,~,~] = DTQP_solve(setup,opts);

%% plots
wcolor = [1 1 1]; % white color
bcolor = [0 0 0]; % black color
fontlabel = 20; % x,y label font size
fontlegend = 12; % x,y legend font size
fonttick = 12; % x,y tick font size

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

%--- state
hf = figure; % create a new figure and save handle
hf.Color = [1 1 1]; % change the figure background color
hf.Position = [200 200 550 400]; % set figure size and position

t = linspace(p.t0,p.tf,10000);

plot(tzo,Xzo); hold on
plot(tef,Xef); hold on
plot(ttr,Xtr); hold on
plot(ths,Xhs); hold on
plot(ths,Xrk4); hold on
plot(tlgl,Xlgl); hold on
plot(t,F(t,p.x0),'k','linewidth',2); hold on

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

hl = legend('$\xi_{ZO}$','$\xi_{EF}$','$\xi_{TR}$','$\xi_{HS}$','$\xi_{RK4}$',...
    '$\xi_{PS}$','$\xi_{actual}$','location','Best');
hl.FontSize = fontlegend; % change legend font size
hl.EdgeColor = bcolor; % change the legend border to black (not a dark grey)

%--- errors
hf = figure; % create a new figure and save handle
hf.Color = [1 1 1]; % change the figure background color
hf.Position = [200 200 550 400]; % set figure size and position

Ezo = abs(Xzo - F(tzo,p.x0));
Eef = abs(Xef - F(tef,p.x0));
Etr = abs(Xtr - F(ttr,p.x0));
Ehs = abs(Xhs - F(ths,p.x0));
Erk4 = abs(Xrk4 - F(trk4,p.x0));
Elgl = abs(Xlgl - F(tlgl,p.x0));

semilogy(tzo,Ezo+1e-20); hold on
semilogy(tef,Eef+1e-20); hold on
semilogy(ttr,Etr+1e-20); hold on
semilogy(ths,Ehs+1e-20); hold on
semilogy(trk4,Erk4+1e-20); hold on
semilogy(tlgl,Elgl+1e-20); hold on

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

hl = legend('$E_{ZO}$','$E_{EF}$','$E_{TR}$','$E_{HS}$','$E_{RK4}$','$E_{PS}$');
hl.FontSize = fontlegend; % change legend font size
hl.EdgeColor = bcolor; % change the legend border to black (not a dark grey)