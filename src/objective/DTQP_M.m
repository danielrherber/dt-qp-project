%--------------------------------------------------------------------------
% DTQP_M.m
% Create sequences for Mayer terms
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [I,J,V] = DTQP_M(Mfull,p,opts)

% initialize sequences  
I = []; J = []; V = []; 

for k = 1:length(Mfull)
    % obtain current substructure
    M = Mfull(k);
    
    % check if it is a matrix
    if ~ismatrix(M.matrix), Mt = cell2mat(M.matrix); else Mt = M.matrix; end
    
    % check if both left and right fields are present, assign 0 if not
    if ~isfield(M,'left'), M.left = 0; end
    if ~isfield(M,'right'), M.right = 0; end
    
    % check if both left and right fields are equal to 0
    if (M.left ~= 0), R = p.i{M.left}; else R = 0; end
    if (M.right ~= 0), C = p.i{M.right}; else C = 0; end
    
    % determine locations and matrix values at this points
    for i = 1:length(R) % number of row continuous variables
        for j = 1:length(C) % number of column continuous variables
            
            r = DTQP_getQPIndex(R(i),M.left,0,p); % Hessian row index sequence
            c = DTQP_getQPIndex(C(j),M.right,0,p); % Hessian column index sequence 
            
            I = [I;r]; % main diagonal row
            J = [J;c]; % main diagonal column

            V = [V;squeeze(Mt(i,j,:))]; % combine
            
        end % end C for loop
    end % end R for loop
end % end M for loop

end % end DTQP_M function