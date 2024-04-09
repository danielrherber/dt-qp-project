%--------------------------------------------------------------------------
% EarthLaunch_plot.m
% Plot function for EarthLaunch example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function EarthLaunch_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% scale time
T = P(1)*T;

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

%% altitude vs. range
figure('Color',wcolor); hold on

% plot
plot(Y(:,1)/1000,Y(:,2)/1000,'.','markersize',12,'color',cArray(1,:));

% highlight final point
ha = gca;
plot(Y(end,1)/1000,Y(end,2)/1000,'k.','markersize',20)
line([ha.XLim(1) Y(end,1)/1000],[Y(end,2) Y(end,2)]/1000,...
    'linestyle','--','Color','k','LineWidth',1);
text(ha.XLim(1)+1,Y(end,2)/1000-6,string(Y(end,2)/1000))
line([Y(end,1) Y(end,1)]/1000,[ha.YLim(1) Y(end,2)/1000],...
    'linestyle','--','Color','k','LineWidth',1);
text(Y(end,1)/1000-15,ha.YLim(1)+6,string(Y(end,1)/1000))

% axis
xlabel('range [km]','fontsize',16)
ylabel('altitude [km]','fontsize',16)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-altitude-range'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end