%--------------------------------------------------------------------------
% DTQP_sensmethods.m
% Solve a specified example using multiple methods and values of nt
%--------------------------------------------------------------------------
% NOTE: only works if you use the standardized output format
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP_sensmethods

% example function
f = @BrysonHo166; % <- change this for different examples
% f = @BrysonDenham;

% create nt test array
N = 100;
ntmin = 4;
ntmax = 1000;
narray = GenerateNgrid(N,ntmin,ntmax);

% method lists
DefectArray = {'ZO','EF','TR','HS','HS','RK4','PS','PS'};
QuadArray = {'CEF','CEF','CTR','CTR','CQHS','CQHS','G','CC'};
MeshArray = {'ED','ED','ED','ED','ED','ED','LGL','CGL'};
LegendArray = strcat(DefectArray,'-',QuadArray,'-',MeshArray);

% select the methods to test from lists above
Testarray = 1:8; % test all methods
% Testarray = [5,6]; % only test CQHS methods
% Testarray = [7,8]; % only test PS methods

DefectArray = DefectArray(Testarray);
QuadArray = QuadArray(Testarray);
MeshArray = MeshArray(Testarray);
LegendArray = LegendArray(Testarray);

% maximum N for PS methods
PSntmax = 500;

% calculate the number of outputs
opts = options();
opts.dt.nt = 5; % dummy value
Ot = f([],opts);
nO = length(Ot);

% grid of test cases
[I,J] = ndgrid(1:length(narray),1:length(Testarray));

% number of tests
Nparfor = numel(I);
Ntestarray = length(Testarray);

% initialize outputs with NaN
output = Ot;
poutput = nan(Nparfor,nO);

% random ordering
QuadV = QuadArray(J);
DefectV = DefectArray(J);
MeshV = MeshArray(J);
NV = narray(I);

% parfor idx = 1:Nparfor % parallel
for idx = 1:Nparfor % serial
    % initialize
    opts = options();
    p = [];

    % 
    opts.dt.quadrature = QuadV{idx};
    opts.dt.defects = DefectV{idx};
    opts.dt.mesh = MeshV{idx}; 
    opts.dt.nt = NV(idx); 

    % check if 
    if sum(strcmpi(opts.dt.mesh,{'LGL','CGL'})) && (opts.dt.nt > PSntmax)
        % do nothing, nan already assigned
    else
        try
            % run and get outputs
            O = f(p,opts);

            % assign outputs to the output matrix
            for k = 1:nO
                poutput(idx,k) = O(k).value;
            end

        catch
            % do nothing, nan already assigned
        end
    end
end

% assign results to a structure
for k = 1:nO
    output(k).value = reshape(poutput(:,k),[],Ntestarray);
end

% plot the results
close all
for k = 1:nO
    outputPlot(output(k),narray,LegendArray);
end

end
% plotting function
function outputPlot(output,narray,LegendArray)
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
fontsize = 16;
map = parula(length(LegendArray));
figure('color',[1 1 1])
for k = 1:length(LegendArray)
    % data
    y = output.value(:,k);
    % check if it is a time calculation
    if contains(output.label,'time','IgnoreCase',true)
        % plot
        loglog(narray,y,'color',map(k,:),'linewidth',2); hold on
    else
        % forward-looking moving maximum
        y = movmax(y,[0 3]);
        % plot
        loglog(narray,y,'color',map(k,:),'linewidth',2); hold on
    end
end
% axis
xlim([narray(1) narray(end)])
xlabel('$n_t$','fontsize',fontsize)
ylabel(['output: ',output.label],'fontsize',fontsize)
% legend
legend(LegendArray,'location','best','fontsize',fontsize-8)
end
% create nt test array
function narray = GenerateNgrid(n,Nmin,Nmax)
% logarithmically spaced integers
Narray = round(logspace(log10(Nmin),log10(Nmax),n));
% random to prevent sequentially even or odd grids
Narray = unique(Narray);
for idx = 1:length(Narray)-3
    EO = mod(Narray(idx),2);
    if all(EO==mod(Narray(idx+1:idx+3),2))
        Narray(idx+1) = Narray(idx+1) + (randi([0,1], 1)*2 - 1);
    end
end
% remove values outside the defined boundaries
Narray(Narray<Nmin) = [];
Narray(Narray>Nmax) = [];
% test only the unique values
narray = unique([Nmin,Narray,Nmax]);
end
% default options
function opts = options
opts.general.plotflag = 0;
opts.general.saveflag = 0;
opts.general.displevel = 1;
opts.qp.disp = 'none';
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED'; 
end