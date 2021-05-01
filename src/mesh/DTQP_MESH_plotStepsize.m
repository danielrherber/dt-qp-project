%--------------------------------------------------------------------------
% DTQP_MESH_plotStepsize.m
% plot the stepsizes for the different meshes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP_MESH_plotStepsize(T)

% plot customizations
fonttick = 12; % x,y tick font size
fontlabel = 20; % x,y label font size
bcolor = [0 0 0]; % black color

% change the interpreters to latex
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');

% initialize plot
hf = figure; hold on; hf.Color = 'w';

% colors
c = parula(length(T));

% initialize min and max stepsize
hmaxplot = -inf; hminplot = inf;

% go through each mesh
for k = 1:length(T)

    % extract
    x = T{k};
    y = diff(T{k});
    y = [y;y(end)];

    % update min and max values
    hmaxplot = max(hmaxplot,max(y));
    hminplot = min(hminplot,min(y));

    % copy for fill command
    x = [x(1);repelem(x(2:end),2)];
    y = [repelem(y(1:end-1),2);y(end)];

    % plot the mesh errors
    fill([x;flipud(x)],[log10(y);-15*ones(size(y))],c(k,:),...
        'FaceAlpha',0.90);

end

% axis labels
xlabel('$t$')
ylabel('Stepsize')

% change axis properties
ha = gca;
ha.YLim = log10([hminplot hmaxplot]);
ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
ha.YAxis.Color = bcolor; % change the y axis color to black (not a dark grey)
ha.XAxis.FontSize = fonttick; % change x tick font size
ha.YAxis.FontSize = fonttick; % change y tick font size
ha.XAxis.Label.FontSize = fontlabel; % change x label font size
ha.YAxis.Label.FontSize = fontlabel; % change y label font size
ha.Layer = 'top'; % place the axes on top of the data
ytickformat('$10^{%.3g}$')

end