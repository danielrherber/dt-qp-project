%--------------------------------------------------------------------------
% BrysonHo63_plot.m
% Plot function for BrysonHo63 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function BrysonHo63_plot(T,U,Y,P,F,in,opts,sol)

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
figname = 'figure-state'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
solutionflag = false; %#ok<NASGU>
flag = 'plot-control'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% x-y
figure('Color',wcolor); hold on

% plot x-y (true)
E = sol(2).E;
plot(E(:,1),E(:,2),'linewidth',2,'color','k')

% plot x-y (from DTQP)
plot(Y(:,1),Y(:,2),'linewidth',2,'color','r')

% axis
xlabel('$x$','fontsize',fontsize_)
ylabel('$y$','fontsize',fontsize_)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>
axis equal

% legend
Lv{1} = '$(x/y)^{*}$';
Lv{2} = '$(x/y)^{DT}$';
hL = legend(Lv);
set(hL,'interpreter','latex','location','best',...
    'fontsize',fontsize_-4)

% save
figname = 'figure-xy'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end