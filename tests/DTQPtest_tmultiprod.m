%--------------------------------------------------------------------------
% DTQPtest_tmultiprod.m
% Testing the DTQP_tmultiprod function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

%--------------------------------------------------------------------------
% single constant scalar
A = rand(1);
P = A;
matrices = A;
p.t = (0:4)';

A = DTQP_tmultiprod(matrices,p);

norm(squeeze(A(1,:,:))-P,'inf')

%--------------------------------------------------------------------------
% single constant matrix
A = rand(6);
P = A;
matrices = A;
p.t = (0:4)';

A = DTQP_tmultiprod(matrices,p);

norm(squeeze(A(1,:,:))-P,'inf')

%--------------------------------------------------------------------------
% single constant matrix (cell)
A = rand(6);
P = A;
matrices = {A};
p.t = (0:4)';

A = DTQP_tmultiprod(matrices,p);

norm(squeeze(A(1,:,:))-P,'inf')

%--------------------------------------------------------------------------
% single time-varying matrix 
A = @(t) 2*exp(-t);
matrices = A;
p.t = (0:4)';
P = 2*exp(-p.t);

A = DTQP_tmultiprod(matrices,p);

norm(A-P,'inf')

%--------------------------------------------------------------------------
% single time-varying matrix (cell)
A = {@(t) 2*exp(-t)};
matrices = A;
p.t = (0:4)';
P = 2*exp(-p.t);

A = DTQP_tmultiprod(matrices,p);

norm(A-P,'inf')

%--------------------------------------------------------------------------
% two constant matrices
A1 = rand(2,2);
A2 = rand(2,2);
P = A1*A2;
matrices = {'prod',A1,A2};
p.t = (0:4)';

A = DTQP_tmultiprod(matrices,p);

norm(squeeze(A(1,:,:))-P,'inf')

%--------------------------------------------------------------------------
% two time-varying matrices
clear
A1 = {@(t) sin(t),0;1,0};
A2 = {@(t) exp(-t) + 1,1;0,0};
matrices = {'prod',A1,A2};
p.t = (0:100)';

A = DTQP_tmultiprod(matrices,p);

P = zeros(length(p.t),2,2);
for k = 1:length(p.t)
    t = p.t(k);
    P(k,:,:) = [sin(t),0;1,0]*[exp(-t) + 1,1;0,0];
end

norm(A(:)-P(:),'inf')

%--------------------------------------------------------------------------
% multiple input matrices
clear
A1 = {@(t) sin(t),0;1,0};
A2 = {1,1;0,0};
A3 = rand(2,10);
A4 = ones(10,1);
A5 = -2;
matrices = {'prod',A1,A2,A3,A4,A5};
p.t = (0:100)';

A = DTQP_tmultiprod(matrices,p);