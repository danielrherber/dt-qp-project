%--------------------------------------------------------------------------
% MultiphaseParameter_plot.m
% Plot function for MultiphaseParameter example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function MultiphaseParameter_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% line colors
cArray = parula(size(Y,2));

% plot state
for i = 1:size(Y,2)
    plot(T,Y(:,i),'.','color',cArray(i,:),'markersize',12);
    plot(sol(2).T,sol(2).Y(:,i),'linewidth',2,'color',cArray(i,:));
end

% axis
xlabel('$t$ (s)')
ylabel('$\xi$')

% legend
Lv = {};
for i = 1:size(Y,2)
    Lv{end+1} = ['$\xi^{DT}_{',num2str(i),'}$'];
    Lv{end+1} = ['$\xi^*_{',num2str(i),'}$'];
end
hL = legend(Lv{1:2}); % only the first two entries
set(hL,'interpreter','latex','location','best','fontsize',fontsize-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% error
figure('Color',wcolor); hold on

% plot error
flag = 'plot-error'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-error'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end