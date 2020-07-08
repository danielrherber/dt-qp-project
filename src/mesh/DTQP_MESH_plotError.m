%--------------------------------------------------------------------------
% DTQP_MESH_plotError.m
% plot the discretization errors for the different meshes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP_MESH_plotError(T,E,tol)

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

% min and max error
emaxplot = log10(max(vertcat(E{:})))+1;
eminplot = log10(min(vertcat(E{:})))-1;

% go through each mesh
for k = 1:length(T)
    % extract
    x = T{k};
    y = [E{k};E{k}(end)];

    % copy for fill command
    x = [x(1);repelem(x(2:end),2)];
    y = [repelem(y(1:end-1),2);y(end)];

    % plot the mesh errors
    fill([x;flipud(x)],[log10(y);eminplot*ones(size(y))],c(k,:),...
        'FaceAlpha',0.90);
end

% plot the tolerance line
plot([T{1}(1) T{1}(end)],log10([tol tol]),'k--')

% axis labels
xlabel('$t$')
ylabel('Discretization Error')

% change axis properties
ha = gca;
ha.YLim = [eminplot emaxplot];
ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark grey)
ha.YAxis.Color = bcolor; % change the y axis color to black (not a dark grey)
ha.XAxis.FontSize = fonttick; % change x tick font size
ha.YAxis.FontSize = fonttick; % change y tick font size
ha.XAxis.Label.FontSize = fontlabel; % change x label font size
ha.YAxis.Label.FontSize = fontlabel; % change y label font size
ha.Layer = 'top'; % place the axes on top of the data
ytickformat('$10^{%.3g}$')
end