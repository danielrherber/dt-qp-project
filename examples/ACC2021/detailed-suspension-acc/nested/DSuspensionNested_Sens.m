%--------------------------------------------------------------------------
% DSuspensionNested_Sens.m
% Implementation sensitivity study using the nested coordination strategy
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% list options to test
x1 = [false,true]; % GA flag
x2 = {'forward','central'}; % derivative method
x3 = [200,600,2000]; % nt
x4 = logspace(-4,-12,3); % inner-loop tolerance
x5 = logspace(-1,-5,3); % outer-loop tolerance

% generate all permutations
[X1,X2,X3,X4,X5] = ndgrid(x1,x2,x3,x4,x5);

% number of elements
N = numel(X1);

% initialize
SOL = cell(N,1);

%--- run tests
% start timer
t2 = tic;

% go through each test
for k = 1:N

    % assign options
    p.GAflag = X1(k);
    p.FiniteDifferenceType = X2{k};
    p.nt = X3(k);
    p.InnerLoopTolerance = X4(k);
    p.OuterLoopTolerance = X5(k);

    % try to solve the problem using nested
    try
        sol = DSuspensionNested(p);
    catch
        sol = nan;
    end

    % assign
    SOL{k} = sol;

    % display progress
    disp(strcat(string(k),"/",string(N)))

end

% end timer
toc(t2)