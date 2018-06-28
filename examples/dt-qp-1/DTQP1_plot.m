%--------------------------------------------------------------------------
% DTQP1_plot.m
% Plot function for DTQP1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP1_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract parameter structure
p = in.p;

close all

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

fontsize = 16;

%% state
figure('Color',[1 1 1]);

% plot state
cArray = lines(size(Y,2));
for i = 1:3
    plot(T,Y(:,i),'.-','color',cArray(i,:),'markersize',12,'linewidth',2); hold on
end

plot([T(1) T(end)], [Y(1,2) Y(end,2)],'ok','markersize',8); hold on
plot(T, p.g(T),'--','color',[0 0 0]); hold on
plot(T, P*ones(size(T)),'--','color',[0.5 0.5 0.5]); hold on

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$\xi$','fontsize',fontsize)

% legend
hL = legend('$\xi_1^*$','$\xi_2^*$','$\xi_3^*$',...
    '$\xi_2^*(t_0) - \xi_2^*(t_f) = 0$',...
    '$\xi_2^* \leq g(t)$',...
    '$\max(\xi_3^*) \leq p^*$' );
set(hL,'location','best','fontsize',fontsize-4,'box','on')

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-state'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

%% control
figure('Color',[1 1 1]);

% plot control
cArray = lines(size(U,2));
for i = 1:size(U,2)
    plot(T,U(:,i),'.-','color',cArray(i,:),'markersize',12,'linewidth',2); hold on
end
plot(T, -10*ones(size(T)),'--','color',[0 0 0]); hold on
plot(T, 10*ones(size(T)),'--','color',[0 0 0]); hold on

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$u$','fontsize',fontsize)

% legend
hL = legend('$u_1^*$','$u_2^*$','$|u_2^*|\leq 10$');
set(hL,'location','best','fontsize',fontsize-4,'box','on')

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-control'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

%% integral and mixed state-control path constraint
figure('Color',[1 1 1]);

% plot
cArray = lines(2);
plot(T,Y(:,4),'.-','linewidth',2,'markersize',12,'color',cArray(1,:)); hold on
plot(T,-Y(:,1)+U(:,2)/12,'.-','linewidth',2,'markersize',12,'color',cArray(2,:)); hold on
plot(T(1), Y(1,4),'sk','markersize',8); hold on
plot(T(end), Y(end,4),'sk','markersize',8); hold on
plot(T, 0*ones(size(T)),'--','color',[0 0 0]); hold on

% axis
xlabel('$t$ (s)','fontsize',fontsize)
xlim([in.t0 in.tf])
ylim([-2.1 0.2])

% legend
hL = legend('$\int_0^t \xi_1^*(s) ds$',...
    '$-\xi^*_1+u^*_2/12 \leq 0$',...
    '$\int_0^1 \xi^*_1(t) dt = 0$');
set(hL,'location','best','fontsize',fontsize-4,'box','on')
hL.Position = [0.3905041 0.372839 0.243739 0.188766];

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-other'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

end