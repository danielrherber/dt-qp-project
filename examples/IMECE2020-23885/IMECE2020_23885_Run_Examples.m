%--------------------------------------------------------------------------
% IMECE2020_23885_Run_Examples.m
% Paper coming soon
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% examples to test
f = @ContainerCrane;
% f = @Vanderpol;
% f = @SimpleCoDesignTransfer;
% f = @SimpleSuspensionSimultaneous;

% DTQP options test cases defined below
testcase = "TR";
% testcase = "PS";
% testcase = "TR-Suspension";

% method fragments
Methods = {'qlin','nonlinearprogram'};
Derivatives = {'symbolic','real-forward'};
Defects = {'TR','PS'};
Quads = {'CTR','G'};
Meshes = {'ED','LGL'};

switch testcase
    %----------------------------------------------------------------------
    case "TR"
    % create nt test array
    narray = [20 200 2000]; % TR

    % method lists
    MethodArray = Methods([1 2 2 2 2]);
    OLQArray = {true,true,false,true,false};
    DerivArray = Derivatives([1 1 1 2 2]);
    DefectArray = Defects([1 1 1 1 1]);
    QuadArray = Quads([1 1 1 1 1]);
    MeshArray = Meshes([1 1 1 1 1]);

    % options
    opts.solver.tolerance = 1e-12;
    %----------------------------------------------------------------------
    case "PS"
    % create nt test array
    narray = [10 40]; % PS

    % method lists
    MethodArray = Methods([1 2 2 2 2]);
    OLQArray = {true,true,false,true,false};
    DerivArray = Derivatives([1 1 1 2 2]);
    DefectArray = Defects([2 2 2 2 2]);
    QuadArray = Quads([2 2 2 2 2]);
    MeshArray = Meshes([2 2 2 2 2]);

    % options
    opts.solver.tolerance = 1e-12;
    %----------------------------------------------------------------------
    case "TR-Suspension" % remove qlin test
    % create nt test array
    narray = [200 2000]; % TR

    % method lists
    MethodArray = Methods([2 2 2 2]);
    OLQArray = {true,false,true,false};
    DerivArray = Derivatives([1 1 2 2]);
    DefectArray = Defects([1 1 1 1]);
    QuadArray = Quads([1 1 1 1]);

    % options
    opts.solver.tolerance = 1e-7;
    %----------------------------------------------------------------------
end

% select the methods to test from lists above
Testarray = 1:length(DefectArray);

% maximum N for PS methods
PSntmax = 400;

% calculate the number of outputs
opts = options(opts);
opts.dt.nt = 200; % dummy value
Ot = f([],opts);
nO = length(Ot);

% grid of test cases
[I,J] = ndgrid(1:length(narray),1:length(Testarray));

% transpose so nt are grouped
I = I'; J = J';

% number of tests
Nparfor = numel(I);
Ntestarray = length(Testarray);

% initialize outputs with NaN
output = Ot;
poutput = nan(Nparfor,nO);

% create matrices with all permutations to test
MethodPermutations = MethodArray(J);
OLQPermutations = OLQArray(J);
DerivPermutations = DerivArray(J);
QuadPermutations = QuadArray(J);
DefectPermutations = DefectArray(J);
MeshPermutations = MeshArray(J);
NPermutations = narray(I);

% parfor idx = 1:Nparfor % parallel
for idx = 1:Nparfor % serial

    % initialize options
    opts = options(opts);

    % extract options set above
    opts.dt.quadrature = QuadPermutations{idx};
    opts.dt.defects = DefectPermutations{idx};
    opts.dt.mesh = MeshPermutations{idx};
    opts.dt.nt = NPermutations(idx);
    opts.method.olqflag = OLQPermutations{idx};
    opts.method.derivatives = DerivPermutations{idx};
    opts.method.form = MethodPermutations{idx};

    % check if
    if sum(strcmpi(opts.dt.mesh,{'LGL','CGL'})) && (opts.dt.nt > PSntmax)
        % do nothing, nan already assigned
    else
        try
            % run and get outputs
            O = f([],opts);

            % assign outputs to the output matrix
            for k = 1:nO
                poutput(idx,k) = O(k).value;
            end

        catch
            % do nothing, nan already assigned
        end
    end

    % display result
    str = string();

    % main method
    switch MethodPermutations{idx}
        %------------------------------------------------------------------
        case "nonlinearprogram"
        str(1) = "\textsl{IP-";
        %------------------------------------------------------------------
        case "qlin"
        str(1) = "\textsl{QLIN-";
        %------------------------------------------------------------------
    end

    % derivative method
    switch DerivPermutations{idx}
        %------------------------------------------------------------------
        case "symbolic"
        str(end+1) = "SYM-";
        %------------------------------------------------------------------
        case "real-forward"
        str(end+1) = "FD-";
        %------------------------------------------------------------------
    end

    % OLQ flag
    switch OLQPermutations{idx}
        %------------------------------------------------------------------
        case true
        str(end+1) = "OLQ-";
        %------------------------------------------------------------------
        case false
        str(end+1) = "NL-";
        %------------------------------------------------------------------
    end

    % defect method and number of points
    str(end+1) = DefectPermutations{idx};
    str(end+1) = string(NPermutations(idx));
    str(end+1) = "}";

    % v
    str(end+1) = " & ";
    str(end+1) = string(sprintf('%0.4f',O(1).value));

    % Iter
    str(end+1) = " & ";
    str(end+1) = string(sprintf('%i',O(6).value));

    % Tsym
    str(end+1) = " & ";
    str(end+1) = string(sprintf('%0.2f',O(2).value));

    % Tint
    str(end+1) = " & ";
    str(end+1) = string(sprintf('%0.3f',O(3).value));

    % Topt
    str(end+1) = " & ";
    str(end+1) = string(sprintf('%0.2f',O(4).value));

    % T
    str(end+1) = " & ";
    str(end+1) = string(sprintf('%0.2f',O(5).value));

    % final return
    str(end+1) = " \\ ";

    % display
    disp(strjoin(str,''))

end

% assign results to a structure
for k = 1:nO
    output(k).value = reshape(poutput(:,k),[],Ntestarray);
end

% default options
function opts = options(opts)

opts.general.plotflag = 0;
opts.general.saveflag = 0;
opts.general.displevel = 1;
opts.solver.display = 'none';
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.method.form = 'nonlinearprogram';
opts.method.olqflag = true;
opts.solver.maxiters = 4000;
opts.method.trustregionflag = false;
opts.method.improveguess = false;

end