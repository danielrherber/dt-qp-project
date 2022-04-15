%--------------------------------------------------------------------------
% TransferMinFuel_plot.m
% Plot function for TransferMinFuel example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function TransferMinFuel_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% plot state
solutionflag = false; %#ok<NASGU>
flag = 'plot-state'; DTQP_plotCommon %#ok<NASGU>

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

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% r vs. theta
figure('Color',wcolor); clf;

polarplot(Y(:,2),Y(:,1),'linewidth',2,'color','r'); hold on
polarplot(Y(1,2),Y(1,1),'k.','markersize',24);
polarplot(Y(end,2),Y(end,1),'bo','markersize',12);

ha = gca;
ha.TickLabelInterpreter = 'latex';
ha.ThetaAxis.Label.String = '$\theta$';
ha.RAxis.Label.String = '$r$';
ha.ThetaAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
ha.RAxis.Color = bcolor; % change the y axis color to black (not a dark grey)
ha.ThetaAxis.FontSize = fontsize_-6; % change x tick font size
ha.RAxis.FontSize = fontsize_-6; % change y tick font size
ha.ThetaAxis.Label.FontSize = fontsize_; % change x label font size
ha.RAxis.Label.FontSize = fontsize_; % change y label font size
ha.Box = 'on'; % box on
ha.Layer = 'top'; % place the axes on top of the data
ha.LineWidth = 1; % increase axis line width

end