%--------------------------------------------------------------------------
% LQRInhomogeneous_plot.m
% Plot function for LQRInhomogeneous example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function LQRInhomogeneous_plot(T,U,Y,P,F,in,opts,sol)

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
    plot(sol(2).T,sol(2).U(:,i),'linewidth',2,'color',cArray(i,:));
end

% axis
xlabel('$t$ (s)')
ylabel('$u$')

% legend
Lv = {};
for i = 1:size(U,2)
    Lv{end+1} = ['$u^{DT}_{',num2str(i),'}$'];
    Lv{end+1} = ['$u^*_{',num2str(i),'}$'];
end
hL = legend(Lv{1:2});
set(hL,'interpreter','latex','location','best','fontsize',fontsize_-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% error
figure('Color',wcolor); hold on

% initialize
Hv = []; Lv = {};

% plot state error
ph1 = semilogy(T,max(abs(sol(1).Y-Y),[],2)+1e-20,'linewidth',2);
Lv{1} = ['$\mathrm{log}_{10}\max_k\|\xi^*_{k}(t) - \xi^{DT}_{k}(t) \|$'];

% plot control error
ph1 = semilogy(T,max(abs(sol(1).U-U),[],2)+1e-20,'linewidth',2);
Lv{2} = ['$\mathrm{log}_{10}\max_k\|u^*_{k}(t) - u^{DT}_{k}(t) \|$'];

% axis
ylim([1e-17 1e0])
xlabel('$t$ (s)')
ylabel('$\mathrm{log}_{10}$ error')
ha = gca;
ha.YScale = 'log';

% legend
hL = legend(Hv,Lv);
set(hL,'location','best','fontsize',fontsize_-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-error'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% mesh refinement
if isfield(in,'meshr')

    % step size plot
    DTQP_MESH_plotStepsize(in.meshr.T)

    % mesh error plot
    DTQP_MESH_plotError(in.meshr.T,in.meshr.errors,in.meshr.tol)

end

end