%--------------------------------------------------------------------------
% DSuspension_plot.m
% Plot function for DetailedSuspension example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DSuspension_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% preliminary plot options
fontsize = 16;
wcolor = [1 1 1]; % white color
bcolor = [0 0 0]; % black color
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

%% state
figure('Color',wcolor); hold on

% plot state
solutionflag = false; %#ok<NASGU>
flag = 'plot-state'; DTQP_plotCommon %#ok<NASGU>

% change line styles
ha = gca;
ha.Children(1).LineStyle = '-';
ha.Children(2).LineStyle = '-';
ha.Children(3).LineStyle = '-';
ha.Children(4).LineStyle = '-';
ha.Children(1).LineWidth = 1.5;
ha.Children(2).LineWidth = 1.5;
ha.Children(3).LineWidth = 1.5;
ha.Children(4).LineWidth = 1.5;

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
solutionflag = false; %#ok<NASGU>
flag = 'plot-control'; DTQP_plotCommon %#ok<NASGU>

% change line styles
ha = gca;
ha.Children(1).LineStyle = '-';
ha.Children(1).LineWidth = 1.5;

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% zs - z0
figure('Color',wcolor); hold on

% line colors
cArray = lines(1);

% plot zs - z0
plot(T,Y(:,1)+Y(:,3),'.-','color',cArray(1,:),'markersize',12,'linewidth',1.5);

% axis
xlabel('$t$ (time)','fontsize',fontsize)
ylabel('$y_1+y_3$','fontsize',fontsize)

% save
figname = 'figure-y1+y3'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end