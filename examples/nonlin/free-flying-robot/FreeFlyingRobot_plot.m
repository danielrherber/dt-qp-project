%--------------------------------------------------------------------------
% FreeFlyingRobot_plot.m
% Plot function for FreeFlyingRobot example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function FreeFlyingRobot_plot(T,U,Y,P,F,in,opts,sol)

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

%% control 2
figure('Color',wcolor); hold on

% line colors
cArray = lines(size(Y,2));

% plot control
plot(T,U(:,1)-U(:,2),'.-','color',cArray(1,:),'markersize',12);
plot(T,U(:,3)-U(:,4),'.-','color',cArray(2,:),'markersize',12);

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$T$','fontsize',fontsize)

% legend
Lv = {};
Lv{end+1} = ['$T^{DT}_{',num2str(1),'}$'];
Lv{end+1} = ['$T^{DT}_{',num2str(2),'}$'];
hL = legend(Lv);
set(hL,'interpreter','latex','location','best',...
    'fontsize',fontsize-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control2'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end