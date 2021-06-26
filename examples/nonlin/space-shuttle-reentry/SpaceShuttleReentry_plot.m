%--------------------------------------------------------------------------
% SpaceShuttleReentry_plot.m
% Plot function for Space Shuttle Reentry example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function SpaceShuttleReentry_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract
auxdata = in.auxdata;

% unscale time
T = T*P;

% calculate heating
vI = 3.280839895*Y(:,4)*auxdata.vs;
hI = 3.280839895*(Y(:,1)*auxdata.rads - auxdata.Re);
rho = 0.002378*exp(-hI/23800);
qr = 17700*sqrt(rho).*(0.0001*vI).^3.07;
a = rad2deg(U(:,1));
qa = 1.0672181 + -0.192137741e-1*a + 0.21286289e-3*a.^2 + -0.10117249e-5*a.^3;
q = qa.*qr;

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

%% states (Betts version)
figure('Color',wcolor); hold on

subplot(3,2,1); hold on
plot(T,(Y(:,1)*auxdata.rads - auxdata.Re)/1000,'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Altitude (1000 m)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,2,3); hold on
plot(T,rad2deg(Y(:,2)),'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Longitude (deg)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,2,5); hold on
plot(T,rad2deg(Y(:,3)),'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Latitude (deg)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,2,2); hold on
plot(T,Y(:,4)*auxdata.vs/1000,'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Velocity (1000 m/s)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,2,4); hold on
plot(T,rad2deg(Y(:,5)),'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Flight Path (deg)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,2,6); hold on
plot(T,rad2deg(Y(:,6)),'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Azimuth (deg)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

% save
figname = 'figure-state-betts'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% controls and heating (Betts version)
figure('Color',wcolor); hold on

subplot(3,1,1); hold on
plot(T,rad2deg(U(:,1)),'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Angle of Attack (deg)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,1,2); hold on
plot(T,rad2deg(U(:,2)),'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Bank Angle (deg)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

subplot(3,1,3); hold on
plot(T,q,'-','linewidth',1.5,'color','k')
xlabel('Time (sec)'); title('Heating (BTU/ft/ft/sec)');
flag = 'axis'; DTQP_plotCommon %#ok<NASGU> % configure axis
ha.XAxis.Label.FontSize = 10;

% save
figname = 'figure-control-betts'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end