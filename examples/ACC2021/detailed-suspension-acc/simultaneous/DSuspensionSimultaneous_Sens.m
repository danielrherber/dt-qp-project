%--------------------------------------------------------------------------
% DSuspensionSimultaneous_Sens.m
% Implementation sensitivity study using the simultaneous coordination
% strategy
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clear; close all; clc

% list options to test
x1 = {'symbolic','complex','real-central','real-forward'}; % derivative method
x2 = [200,600,2000]; % nt
x3 = logspace(-3,-7,3); % optimality tolerance
x4 = logspace(-4,-12,3); % feasibility tolerance

% generate all permutations
[X1,X2,X3,X4] = ndgrid(x1,x2,x3,x4);

% number of elements
N = numel(X1);

% initialize
SOL = cell(N,1);

%--- run tests
% start timer
t2 = tic;

% go through each test
for k = 1:N

    % try to solve the problem using nested
    try
        sol = DSuspensionSimultaneous(X1{k},X2(k),X3(k),X4(k));
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