%--------------------------------------------------------------------------
% DTQP1_plot.m
% Plot function for DTQP1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP1_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract parameter structure
p = in.p;

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% line colors
cArray = lines(size(Y,2));

% plot state
for i = 1:3
    plot(T,Y(:,i),'.-','color',cArray(i,:),'markersize',12,'linewidth',2);
end

plot([T(1) T(end)], [Y(1,2) Y(end,2)],'ok','markersize',8);
plot(T, p.g(T),'--','color',[0 0 0]);
plot(T, P*ones(size(T)),'--','color',[0.5 0.5 0.5]);

% axis
xlabel('$t$ (s)')
ylabel('$\xi$')

% legend
hL = legend('$\xi_1^*$','$\xi_2^*$','$\xi_3^*$',...
    '$\xi_2^*(t_0) - \xi_2^*(t_f) = 0$',...
    '$\xi_2^* \leq g(t)$',...
    '$\max(\xi_3^*) \leq p^*$' );
set(hL,'location','best','fontsize',fontsize-4,'box','on')

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% line colors
cArray = lines(size(Y,2));

% plot control
for i = 1:size(U,2)
    plot(T,U(:,i),'.-','color',cArray(i,:),'markersize',12,'linewidth',2);
end

plot(T, -10*ones(size(T)),'--','color',[0 0 0]);
plot(T, 10*ones(size(T)),'--','color',[0 0 0]);

% axis
xlabel('$t$ (s)')
ylabel('$u$')

% legend
hL = legend('$u_1^*$','$u_2^*$','$|u_2^*|\leq 10$');
set(hL,'location','best','fontsize',fontsize-4,'box','on')

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% integral and mixed state-control path constraint
figure('Color',wcolor); hold on

% plot
cArray = lines(2);
plot(T,Y(:,4),'.-','linewidth',2,'markersize',12,'color',cArray(1,:));
plot(T,-Y(:,1)+U(:,2)/12,'.-','linewidth',2,'markersize',12,'color',cArray(2,:));
plot(T(1), Y(1,4),'sk','markersize',8);
plot(T(end), Y(end,4),'sk','markersize',8);
plot(T, 0*ones(size(T)),'--','color',[0 0 0]);

% axis
xlabel('$t$ (s)')
xlim([in.t0 in.tf])
ylim([-2.1 0.2])

% legend
hL = legend('$\int_0^t \xi_1^*(s) ds$',...
    '$-\xi^*_1+u^*_2/12 \leq 0$',...
    '$\int_0^1 \xi^*_1(t) dt = 0$');
set(hL,'location','best','fontsize',fontsize-4,'box','on')
hL.Position = [0.3905041 0.372839 0.243739 0.188766];

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-other'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end