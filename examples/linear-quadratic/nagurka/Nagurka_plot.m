%--------------------------------------------------------------------------
% Nagurka_plot.m
% Plot function for Nagurka example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Nagurka_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% line colors
cArray = parula(size(Y,2));

% plot state
ny = size(Y,2);
for i = [1,ny,2:ny-1]
    plot(T,Y(:,i),'.','color',cArray(i,:),'markersize',12);
end

% axis
xlabel('$t$ (s)')
ylabel('$\xi$')

% legend
Lv = {};
for i = [1,ny]
    Lv{end+1} = ['$\xi^{DT}_{',num2str(i),'}$'];
end
hL = legend(Lv{:}); % only the first two entries
set(hL,'interpreter','latex','location','best','fontsize',fontsize_-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% line colors
cArray = parula(size(U,2));

% plot control
for i = 1:size(U,2)
    plot(T,U(:,i),'.','color',cArray(i,:),'markersize',12);
end

% axis
xlabel('$t$ (s)')
ylabel('$u$')

% legend
Lv = {};
for i = 1:size(U,2)
    Lv{end+1} = ['$u^{DT}_{',num2str(i),'}$'];
end
hL = legend(Lv{1:2});
set(hL,'interpreter','latex','location','best','fontsize',fontsize_-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end