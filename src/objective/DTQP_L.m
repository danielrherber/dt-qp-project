%--------------------------------------------------------------------------
% DTQP_L.m
% Create sequences for Lagrange terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,J,V] = DTQP_L(Lfull,p,opts)

% number of time points
nt = p.nt;

% do we need off diagonal terms?
OffFlag = any(strcmpi(opts.dt.quadrature,'CQHS'));

% initialize sequences  
I = []; J = []; W = []; H = []; Q = []; Qmid = [];

for k = 1:length(Lfull)
    % obtain current substructure
    L = Lfull(k);

	% find time dependent matrix
	Lt = DTQP_tmultiprod(L.matrix,p);

    % find time dependent matrix at time grid midpoints
    if OffFlag
        Lmid = DTQP_tmultiprod(L.matrix,p,p.tm); 
    end

    % check if both left and right fields are present, assign 0 if not
    if ~isfield(L,'left'), L.left = 0; end
    if ~isfield(L,'right'), L.right = 0; end

    % check if both left and right fields are equal to 0
    if (L.left ~= 0), R = p.i{L.left}; else R = 0; end % rows (continuous)
    if (L.right ~= 0), C = p.i{L.right}; else C = 0; end % columns (continuous)

    % determine locations and matrix values at this points
    for i = 1:length(R) % number of row continuous variables
        for j = 1:length(C) % number of column continuous variables

            r = DTQP_getQPIndex(R(i),L.left,1,p); % Hessian row index sequence
            c = DTQP_getQPIndex(C(j),L.right,1,p); % Hessian column index sequence 

            I = [I;r]; % main diagonal row
            J = [J;c]; % main diagonal column

            % combine integration weights
            if any(strcmpi(opts.dt.quadrature,{'G','CC'}))
                W = [W; p.w];
            else
                H = [H; p.h; 0];
            end

            % combine values evaluated on the time grid
            Q = [Q;Lt(:,i,j)];   

            % combine values evaluated on the intermediate time grid
            if OffFlag
                Lmidv = Lmid(:,i,j);
                Qmid = [Qmid; Lmidv; 0];
            end

        end % end C for loop
    end % end R for loop
end % end L for loop

% calculate off diagonal sequences
if OffFlag
    IU = I; IU(nt:nt:end) = [];
    IL = I; IL(1:nt:end) = [];
    JU = J; JU(1:nt:end) = [];
    JL = J; JL(nt:nt:end) = [];
    Hoff = H; Hoff(nt:nt:end) = [];
    Qoff = Q; Qoff(nt:nt:end) = [];
end

% begin method specific
switch upper(opts.dt.quadrature)
    case 'CEF'
        V = H.*Q;
    case 'CTR'
        V = ( (H.*Q) + (circshift(H,[1,1]).*Q) )/2;
    case 'CQHS'
        V = ( (H.*Q) + (circshift(H,[1,1]).*Q) + (H.*Qmid) + circshift(H.*Qmid,[1,1]) )/6;
        Voff = Hoff.*Qoff/6;
    case 'G'
        V = (p.tf - p.t0)/2*W.*Q;
    case 'CC'
        V = (p.tf - p.t0)/2*W.*Q;
end
% end method specific

% concatenate for final outputs
if OffFlag
    I = [I;IU;IL];
    J = [J;JU;JL];
    V = [V;Voff;Voff];
end

end % end DTQP_L function