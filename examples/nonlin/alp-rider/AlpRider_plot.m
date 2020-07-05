%--------------------------------------------------------------------------
% AlpRider_plot.m
% Plot function for AlpRider example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor : Athul K. Sundarrajan (AthulKrishnaSundarrajan on Github)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function AlpRider_plot(T,U,Y,P,F,in,opts,sol)

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

%% path constraint
figure('Color',wcolor); hold on

% line colors
cArray = lines(size(Y,2));

% plot path constraint
T2 = linspace(T(1),T(end),1e4)';
plot(T2,3*exp(-12*(T2-3).^2) + 3*exp(-10*(T2-6).^2) + 3*exp(-6*(T2-10).^2) + 8*exp(-4*(T2-15).^2) + 0.01,'.','color',cArray(1,:),'markersize',12);
plot(T,sum(Y.^2,2),'.','color',cArray(2,:),'markersize',12);

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('peak function','fontsize',fontsize)

% legend
Lv = {};
Lv{end+1} = ['peak function'];
Lv{end+1} = ['$(\xi_1^2 + \xi_2^2 + \xi_3^2 + \xi_4^2)^{DT}$'];
hL = legend(Lv);
set(hL,'interpreter','latex','location','best',...
    'fontsize',fontsize-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-path'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end