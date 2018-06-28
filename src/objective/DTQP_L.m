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
function [I,J,V] = DTQP_L(Lfull,in,opts)

% extract some of the variables
nt = in.nt; nu = in.nu; ny = in.ny; t = in.t; tm = in.tm; h = in.h;
w = in.w; p = in.p;  ini = in.i; quadrature = opts.dt.quadrature;

% check if we need off-diagonal terms
OffFlag = any(strcmpi(quadrature,'CQHS'));

% initialize storage arrays
Isav = {}; Jsav = {}; HWsav = {}; Qsav = {}; Qmsav = {};
IUsav = {}; ILsav = {}; JUsav = {}; JLsav = {}; Hoffsav = {}; Qoffsav = {};

% go through each Lagrange term
for k = 1:length(Lfull)
    % obtain current substructure
    Lleft = Lfull(k).left;
    Lright = Lfull(k).right;
    Lmatrix = Lfull(k).matrix;

	% find time dependent matrix
	Lt = DTQP_tmultiprod(Lmatrix,p,t);

    % find time dependent matrix at time grid midpoints
    if OffFlag
        Lm = DTQP_tmultiprod(Lmatrix,p,tm); 
    end

    % check if both left and right fields are equal to 0
    if (Lleft ~= 0), R = ini{Lleft}; else, R = 0; end % rows (continuous)
    if (Lright ~= 0), C = ini{Lright}; else, C = 0; end % columns (continuous)

    % determine locations and matrix values at this points
    for i = 1:length(R) % number of row continuous variables
        for j = 1:length(C) % number of column continuous variables
            % get current matrix values
            Lv = Lt(:,i,j);

            % check if this entry is always 0
            if any(Lv)
                r = DTQP_getQPIndex(R(i),Lleft,1,nt,nu,ny); % Hessian row index sequence
                c = DTQP_getQPIndex(C(j),Lright,1,nt,nu,ny); % Hessian column index sequence 

                Isav{end+1} = r; % main diagonal rows
                Jsav{end+1} = c; % main diagonal columns

                % integration weights
                if any(strcmpi(quadrature,{'G','CC'}))
                    HWsav{end+1} = w;
                else
                    HWsav{end+1} = [h; 0];
                end

                % main diagonal matrix values
                Qsav{end+1} = Lv;   

                % off-diagonal sequences
                if OffFlag
                    Lmv = Lm(:,i,j);
                    Qmsav{end+1} = [Lmv; 0];
                    IUsav{end+1} = r(1:end-1);
                    ILsav{end+1} = r(2:end);
                    JUsav{end+1} = c(2:end);
                    JLsav{end+1} = c(1:end-1);
                    Hoffsav{end+1} = h;
                    Qoffsav{end+1} = Lv(1:end-1);
                end

            end

        end % end C for loop
    end % end R for loop
end % end L for loop

% combine sequences
I = vertcat(Isav{:});
J = vertcat(Jsav{:});
HW = vertcat(HWsav{:});
Q = vertcat(Qsav{:});

% combine off-diagonal sequences
if OffFlag
    Qm = vertcat(Qmsav{:});
    IU = vertcat(IUsav{:});
    IL = vertcat(ILsav{:});
    JU = vertcat(JUsav{:});
    JL = vertcat(JLsav{:});
    Hoff = vertcat(Hoffsav{:});
    Qoff = vertcat(Qoffsav{:});
end

% begin method specific
switch upper(opts.dt.quadrature)
    case 'CEF'
        V = HW.*Q;
    case 'CTR'
        V = ( (HW.*Q) + (circshift(HW,[1,1]).*Q) )/2;
    case 'CQHS'
        V = ( (HW.*Q) + (circshift(HW,[1,1]).*Q) + (HW.*Qm) + circshift(HW.*Qm,[1,1]) )/6;
        Voff = Hoff.*Qoff/6;
    case 'G'
        V = (in.tf - in.t0)/2*HW.*Q;
    case 'CC'
        V = (in.tf - in.t0)/2*HW.*Q;
end
% end method specific

% combine
if OffFlag
    I = vertcat(I,IU,IL);
    J = vertcat(J,JU,JL);
    V = vertcat(V,Voff,Voff);
end

end % end DTQP_L function