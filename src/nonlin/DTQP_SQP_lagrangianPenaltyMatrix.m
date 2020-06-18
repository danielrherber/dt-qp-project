%--------------------------------------------------------------------------
% DTQP_SQP_lagrangianPenaltyMatrix.m
% Create sequences for Hessian matrix for the Lagrangian penalty matrix in
% the SQP method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,J,V] = DTQP_SQP_lagrangianPenaltyMatrix(Lmatrix,in,opts)

% extract some of the variables
nu = in.nu; ny = in.ny; np = in.np; ini = in.i;
t = in.t; h = in.h; nt = in.nt; p = in.p;
lambda = opts.lambda;

% extract and reshape multipliers
lambda = lambda.eqlin(in.multipliers.defects);

% augment step size and multiplier vectors with initial zeros
h = [0;h/2];
lambda = vertcat(zeros(1,ny),lambda);

% initialize row and column indices
LR = repelem([1 2 3],[nu ny np]);
R = horzcat(ini{1:3});
C = R;

% initialize storage arrays
Isav = {}; Jsav = {}; Hsav = {}; Lsav = {};  Qsav = {};

% go through each Lagrange term
for k = 1:length(Lmatrix)

    % find time dependent matrix
    Lt = DTQP_tmultiprod(Lmatrix{k},p,t);

    % determine locations and matrix values at this points
    for i = 1:length(R) % number of row continuous variables
        for j = 1:length(C) % number of column continuous variables

            % get current matrix values
            Lv = Lt(:,i,j);

            % check if this entry is always 0
            if any(Lv)
                r = DTQP_getQPIndex(R(i),LR(i),1,nt,nu,ny); % Hessian row index sequence
                c = DTQP_getQPIndex(C(j),LR(j),1,nt,nu,ny); % Hessian column index sequence

                Isav{end+1} = r; % main diagonal rows
                Jsav{end+1} = c; % main diagonal columns

                % step size
                Hsav{end+1} = h;

                % Lagrange multipliers
                Lsav{end+1} = lambda(:,k);

                % main diagonal values
                Qsav{end+1} = -Lv;

            end

        end % end C for loop
    end % end R for loop
end % end Lmatrix for loop

% combine sequences
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
H = vertcat(Hsav{:});
L = vertcat(Lsav{:});
Q = vertcat(Qsav{:});

% compute main diagonal values
LH = L.*H;
V = ((LH + circshift(LH,[-1,0]))).*Q;

end